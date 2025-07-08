/*
#########################################################################################################
# mr_classes_extractor.c
#
# Finds the n values that generate a specific mr value in Collatz sequences using tuple-based transform 
#
# The complete article can be found on https://doi.org/10.5281/zenodo.15546925
#
#########################################################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <omp.h>

/* ============================================================================
 * CONSTANTS AND CONFIGURATION
 * ============================================================================
 */

#define MAX_SEQUENCE_LENGTH 50000
#define INITIAL_M_CAPACITY 100
#define HASH_TABLE_SIZE 8192
#define PROGRESS_UPDATE_INTERVAL 3.0
#define PROGRESS_CHECK_FREQUENCY 100000
#define MIN_CHUNK_SIZE 100
#define MAX_CHUNK_SIZE 10000
#define CHUNK_MULTIPLIER 10
#define INITIAL_RESULTS_CAPACITY 1000
#define RESULTS_FLUSH_INTERVAL 10000

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================
 */

// Hash table node for fast m value lookup
typedef struct HashNode {
    uint64_t value;
    struct HashNode* next;
} HashNode;

// Container for m values with hash table optimization
typedef struct {
    HashNode* buckets[HASH_TABLE_SIZE];
    uint64_t* values;
    int count;
    int capacity;
} MValues;

// Structure to hold detailed taxonomy analysis results
typedef struct {
    char taxonomy;
    int first_mr_pos;
    int second_mr_pos;
    uint64_t max_m_value;
    int max_m_position;
    uint64_t max_collatz_value;
    int max_collatz_position;
    bool found_repetition;
} TaxonomyAnalysis;

// Thread-safe set for n values that generate target mr with taxonomy and detailed info
typedef struct {
    uint64_t* n_values;
    char** taxonomies;
    int* first_mr_pos;
    int* second_mr_pos;
    uint64_t* max_m_value; 
    int* max_m_position; 
    uint64_t* max_collatz_value;
    int* max_collatz_position;
    int count;
    int capacity;
    omp_lock_t lock;
} TargetResults;

// Progress tracking structure
typedef struct {
    uint64_t processed;
    uint64_t found_count;
    uint64_t last_n_found;
    double last_update_time;
    omp_lock_t lock;
} ProgressTracker;

// Main search context containing all state
typedef struct {
    uint64_t max_n;
    uint64_t target_mr;
    int exponent;
    TargetResults* results;
    ProgressTracker* progress;
    double start_time;
    FILE* output_file;
    omp_lock_t file_lock;
} SearchContext;

// Scheduling configuration
typedef enum {
    SCHEDULE_STATIC = 0,
    SCHEDULE_GUIDED = 1,
    SCHEDULE_DYNAMIC = 2
} SchedulingType;

typedef struct {
    SchedulingType type;
    uint64_t chunk_size;
    int num_threads;
} SchedulingConfig;

// Estructura para agrupar resultados por taxonomía
typedef struct {
    char taxonomy[8];
    uint64_t* n_values;
    int count;
    int capacity;
} TaxonomyGroup;

/* ============================================================================
 * FORWARD DECLARATIONS
 * ============================================================================
 */

// TaxonomyGroup function prototypes
static TaxonomyGroup* create_taxonomy_group(const char* taxonomy);
static void add_n_to_group(TaxonomyGroup* group, uint64_t n);
static void destroy_taxonomy_group(TaxonomyGroup* group);
static void sort_n_values_in_group(TaxonomyGroup* group);

// File writing function prototypes
static void write_final_results_to_disk(SearchContext* ctx);

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================
 */

// Hash function for fast lookup
static inline uint64_t hash_function(uint64_t value) {
    return value & (HASH_TABLE_SIZE - 1);
}

// Calculate m parameter for a Collatz value
static inline uint64_t calculate_m(uint64_t c) {
    uint64_t p = (c & 1) ? 1 : 2;
    return (c - p) >> 1;
}

// Apply Collatz transformation with overflow check
static inline bool apply_collatz_transform(uint64_t* n) {
    if (*n & 1) {  // Odd number
        // Check for overflow in 3*n + 1
        if (*n > (UINT64_MAX - 1) / 3) {
            return false;
        }
        uint64_t temp = 3 * (*n);
        if (temp > UINT64_MAX - 1) {
            return false;
        }
        *n = temp + 1;
    } else {  // Even number
        *n = *n >> 1;
    }
    return true;
}

// Safe memory allocation with error checking
static void* safe_malloc(size_t size, const char* context) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Error: Memory allocation failed for %s\n", context);
        exit(1);
    }
    return ptr;
}

// Safe memory reallocation with error checking
static void* safe_realloc(void* ptr, size_t size, const char* context) {
    void* new_ptr = realloc(ptr, size);
    if (!new_ptr) {
        fprintf(stderr, "Error: Memory reallocation failed for %s\n", context);
        exit(1);
    }
    return new_ptr;
}

/* ============================================================================
 * TAXONOMY GROUP FUNCTIONS
 * ============================================================================
 */

// Helper function to create a taxonomy group
static TaxonomyGroup* create_taxonomy_group(const char* taxonomy) {
    TaxonomyGroup* group = safe_malloc(sizeof(TaxonomyGroup), "taxonomy group");
    strcpy(group->taxonomy, taxonomy);
    group->capacity = 100;
    group->n_values = safe_malloc(group->capacity * sizeof(uint64_t), "taxonomy group n_values");
    group->count = 0;
    return group;
}

// Helper function to add an n value to a group
static void add_n_to_group(TaxonomyGroup* group, uint64_t n) {
    if (group->count >= group->capacity) {
        group->capacity *= 2;
        group->n_values = safe_realloc(group->n_values, group->capacity * sizeof(uint64_t), "taxonomy group expansion");
    }
    group->n_values[group->count++] = n;
}

// Helper function to free memory of a group
static void destroy_taxonomy_group(TaxonomyGroup* group) {
    if (group) {
        free(group->n_values);
        free(group);
    }
}

// Helper function to sort n values within a group
static void sort_n_values_in_group(TaxonomyGroup* group) {
    for (int i = 0; i < group->count - 1; i++) {
        for (int j = 0; j < group->count - 1 - i; j++) {
            if (group->n_values[j] > group->n_values[j + 1]) {
                uint64_t temp = group->n_values[j];
                group->n_values[j] = group->n_values[j + 1];
                group->n_values[j + 1] = temp;
            }
        }
    }
}

// Write accumulated results to file and clear memory (only for very large datasets)
static void flush_results_to_disk(SearchContext* ctx) {
    if (ctx->results->count == 0) return;
    
    // Don't flush if we have small to medium number of results - keep in memory for detailed display
    if (ctx->progress->found_count < 50000) {
        return;  // Keep data in memory for display
    }
    
    // For very large datasets, use the final write function and clear memory
    write_final_results_to_disk(ctx);
    
    // Clear results from memory but keep the arrays allocated
    ctx->results->count = 0;
}

// Initialize incremental file writing
static bool init_incremental_output(SearchContext* ctx, const char* filename) {
    ctx->output_file = fopen(filename, "w");
    if (!ctx->output_file) {
        printf("Error: Could not create output file %s\n", filename);
        return false;
    }
    
    // Write CSV header
    fprintf(ctx->output_file, "mr,taxonomy,n_values\n");
    fflush(ctx->output_file);
    
    omp_init_lock(&ctx->file_lock);
    return true;
}

// Close incremental file writing
static void close_incremental_output(SearchContext* ctx) {
    if (ctx->output_file) {
        fclose(ctx->output_file);
        ctx->output_file = NULL;
    }
    omp_destroy_lock(&ctx->file_lock);
}

// Write accumulated results to file (forced write for final output)
static void write_final_results_to_disk(SearchContext* ctx) {
    if (ctx->results->count == 0) return;
    
    omp_set_lock(&ctx->file_lock);
    
    // Group results by taxonomy
    TaxonomyGroup* groups[4] = {NULL}; // A, B, C, Simple
    const char* taxonomy_names[] = {"A", "B", "C", "Simple"};
    
    // Group current results
    for (int i = 0; i < ctx->results->count; i++) {
        const char* taxonomy = ctx->results->taxonomies[i];
        uint64_t n_value = ctx->results->n_values[i];
        
        int group_index = -1;
        for (int j = 0; j < 4; j++) {
            if (strcmp(taxonomy, taxonomy_names[j]) == 0) {
                group_index = j;
                break;
            }
        }
        
        if (group_index != -1) {
            if (groups[group_index] == NULL) {
                groups[group_index] = create_taxonomy_group(taxonomy_names[group_index]);
            }
            add_n_to_group(groups[group_index], n_value);
        }
    }
    
    // Write groups to file
    for (int i = 0; i < 4; i++) {
        if (groups[i] != NULL && groups[i]->count > 0) {
            sort_n_values_in_group(groups[i]);
            
            fprintf(ctx->output_file, "%lu,%s,\"", ctx->target_mr, groups[i]->taxonomy);
            for (int j = 0; j < groups[i]->count; j++) {
                fprintf(ctx->output_file, "%lu", groups[i]->n_values[j]);
                if (j < groups[i]->count - 1) {
                    fprintf(ctx->output_file, ",");
                }
            }
            fprintf(ctx->output_file, "\"\n");
        }
    }
    
    fflush(ctx->output_file);
    
    // Clean up groups
    for (int i = 0; i < 4; i++) {
        if (groups[i] != NULL) {
            destroy_taxonomy_group(groups[i]);
        }
    }
    
    omp_unset_lock(&ctx->file_lock);
}

/* ============================================================================
 * PROGRESS TRACKER IMPLEMENTATION
 * ============================================================================
 */

static ProgressTracker* create_progress_tracker(void) {
    ProgressTracker* tracker = safe_malloc(sizeof(ProgressTracker), "progress tracker");
    tracker->processed = 0;
    tracker->found_count = 0;
    tracker->last_n_found = 0;
    tracker->last_update_time = 0.0;
    omp_init_lock(&tracker->lock);
    return tracker;
}

static void destroy_progress_tracker(ProgressTracker* tracker) {
    if (tracker) {
        omp_destroy_lock(&tracker->lock);
        free(tracker);
    }
}

static void update_last_n_found(ProgressTracker* tracker, uint64_t n) {
    omp_set_lock(&tracker->lock);
    tracker->last_n_found = n;
    omp_unset_lock(&tracker->lock);
}

static void update_progress_if_needed(const SearchContext* ctx) {
    ProgressTracker* tracker = ctx->progress;
    omp_set_lock(&tracker->lock);
    
    double current_time = omp_get_wtime();
    
    if (current_time - tracker->last_update_time >= PROGRESS_UPDATE_INTERVAL) {
        tracker->last_update_time = current_time;
        
        double elapsed = current_time - ctx->start_time;
        double rate = tracker->processed / elapsed;
        double eta = (ctx->max_n - tracker->processed) / rate;
        double progress_percent = (double)tracker->processed / ctx->max_n * 100.0;
        
        printf("Progress: (%.1f%%) | Processed: %lu | Found n values with mr=%lu: %d | %.1f nums/sec | ETA: %.1f min\n",
               progress_percent, tracker->processed, ctx->target_mr, ctx->results->count, rate, eta/60.0);
        fflush(stdout);
    }
    
    omp_unset_lock(&tracker->lock);
}

static void increment_progress_counters(ProgressTracker* tracker, bool target_found) {
    #pragma omp atomic
    tracker->processed++;
    
    if (target_found) {
        #pragma omp atomic
        tracker->found_count++;
    }
}

/* ============================================================================
 * M VALUES CONTAINER IMPLEMENTATION
 * ============================================================================
 */

static void init_m_values(MValues* mv) {
    mv->capacity = INITIAL_M_CAPACITY;
    mv->values = safe_malloc(mv->capacity * sizeof(uint64_t), "m_values array");
    mv->count = 0;
    
    // Initialize hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        mv->buckets[i] = NULL;
    }
}

static void destroy_m_values(MValues* mv) {
    if (!mv) return;
    
    // Clean up hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        HashNode* current = mv->buckets[i];
        while (current) {
            HashNode* next = current->next;
            free(current);
            current = next;
        }
        mv->buckets[i] = NULL;
    }
    
    // Clean up values array
    free(mv->values);
    mv->values = NULL;
    mv->count = 0;
    mv->capacity = 0;
}

static bool is_m_repeated(const MValues* mv, uint64_t m) {
    uint64_t hash = hash_function(m);
    HashNode* current = mv->buckets[hash];
    
    while (current) {
        if (current->value == m) {
            return true;
        }
        current = current->next;
    }
    return false;
}

static void add_m_value(MValues* mv, uint64_t m) {
    // Expand array if necessary
    if (mv->count >= mv->capacity) {
        int new_capacity = mv->capacity * 2;
        mv->values = safe_realloc(mv->values, new_capacity * sizeof(uint64_t), "m_values expansion");
        mv->capacity = new_capacity;
    }
    
    mv->values[mv->count++] = m;
    
    // Add to hash table
    uint64_t hash = hash_function(m);
    HashNode* new_node = safe_malloc(sizeof(HashNode), "hash node");
    new_node->value = m;
    new_node->next = mv->buckets[hash];
    mv->buckets[hash] = new_node;
}

/* ============================================================================
 * TAXONOMY CLASSIFICATION FUNCTIONS
 * ============================================================================
 */

// Classify sequence taxonomy and return detailed analysis
static TaxonomyAnalysis analyze_sequence_taxonomy(uint64_t n_start, uint64_t target_mr) {
    TaxonomyAnalysis result = {0};
    result.taxonomy = 'S';
    result.first_mr_pos = -1;
    result.second_mr_pos = -1;
    result.max_m_value = 0;
    result.max_m_position = -1;
    result.max_collatz_value = 0;
    result.max_collatz_position = -1;
    result.found_repetition = false;
    
    uint64_t n = n_start;
    
    uint64_t* collatz_sequence = safe_malloc(MAX_SEQUENCE_LENGTH * sizeof(uint64_t), "collatz sequence");
    uint64_t* m_sequence = safe_malloc(MAX_SEQUENCE_LENGTH * sizeof(uint64_t), "m sequence");
    int sequence_length = 0;
    int m_sequence_length = 0;
    
    // Build the sequences completely first
    collatz_sequence[sequence_length++] = n;
    
    // Generate the complete m-sequence and look for target_mr repetitions
    for (int step = 0; step < MAX_SEQUENCE_LENGTH - 1 && n != 1; step++) {
        uint64_t m = calculate_m(n);
        m_sequence[m_sequence_length++] = m;
        
        // Apply Collatz transformation
        if (!apply_collatz_transform(&n)) {
            break;
        }
        
        collatz_sequence[sequence_length++] = n;
        
        // Safety check
        if (step > 10000 && n > n_start * 1000) {
            break;
        }
    }
    
    // Now search for target_mr repetitions in the complete m_sequence
    for (int i = 0; i < m_sequence_length; i++) {
        if (m_sequence[i] == target_mr) {
            if (result.first_mr_pos == -1) {
                result.first_mr_pos = i + 1;
            } else {
                // Found second occurrence
                result.second_mr_pos = i + 1;
                result.found_repetition = true;
                break;
            }
        }
    }
    
    // Find maximum value and its position in the Collatz sequence
    if (sequence_length > 0) {
        result.max_collatz_value = collatz_sequence[0];
        result.max_collatz_position = 1;
        
        for (int i = 1; i < sequence_length; i++) {
            if (collatz_sequence[i] > result.max_collatz_value) {
                result.max_collatz_value = collatz_sequence[i];
                result.max_collatz_position = i + 1;
            }
        }
    }
    
    // Find maximum value and its position in the m-parameter sequence
    if (m_sequence_length > 0) {
        result.max_m_value = m_sequence[0];
        result.max_m_position = 1;  // 1-indexed
        
        for (int i = 1; i < m_sequence_length; i++) {
            if (m_sequence[i] > result.max_m_value) {
                result.max_m_value = m_sequence[i];
                result.max_m_position = i + 1;
            }
        }
    }
    
    // Classify taxonomy based on relative positions of MAXIMUM M vs mr positions
    // All positions are now consistently in m_sequence indexing (1-indexed)
    if (result.found_repetition && result.first_mr_pos != -1 && result.second_mr_pos != -1) {
        if (result.max_m_position < result.first_mr_pos) {
            result.taxonomy = 'A';
        } else if (result.max_m_position >= result.first_mr_pos && result.max_m_position < result.second_mr_pos) {
            result.taxonomy = 'B';
        } else {
            result.taxonomy = 'C';
        }
    }
    
    free(collatz_sequence);
    free(m_sequence);
    return result;
}

// Simple wrapper for backward compatibility
static char classify_sequence_taxonomy(uint64_t n_start, uint64_t target_mr) {
    TaxonomyAnalysis analysis = analyze_sequence_taxonomy(n_start, target_mr);
    return analysis.taxonomy;
}

/* ============================================================================
 * TARGET RESULTS IMPLEMENTATION
 * ============================================================================
 */

static TargetResults* create_target_results(void) {
    TargetResults* results = safe_malloc(sizeof(TargetResults), "target results");
    results->capacity = INITIAL_RESULTS_CAPACITY;
    results->n_values = safe_malloc(results->capacity * sizeof(uint64_t), "n values array");
    results->taxonomies = safe_malloc(results->capacity * sizeof(char*), "taxonomies array");
    results->first_mr_pos = safe_malloc(results->capacity * sizeof(int), "first mr pos array");
    results->second_mr_pos = safe_malloc(results->capacity * sizeof(int), "second mr pos array");
    results->max_m_value = safe_malloc(results->capacity * sizeof(uint64_t), "max m value array");
    results->max_m_position = safe_malloc(results->capacity * sizeof(int), "max m position array");
    results->max_collatz_value = safe_malloc(results->capacity * sizeof(uint64_t), "max collatz value array");
    results->max_collatz_position = safe_malloc(results->capacity * sizeof(int), "max collatz position array");
    
    // Initialize taxonomy strings
    for (int i = 0; i < results->capacity; i++) {
        results->taxonomies[i] = safe_malloc(8 * sizeof(char), "taxonomy string");
    }
    
    results->count = 0;
    omp_init_lock(&results->lock);
    return results;
}

static void destroy_target_results(TargetResults* results) {
    if (results) {
        free(results->n_values);
        free(results->first_mr_pos);
        free(results->second_mr_pos);
        free(results->max_m_value);
        free(results->max_m_position);
        free(results->max_collatz_value);
        free(results->max_collatz_position);
        
        // Free taxonomy strings
        for (int i = 0; i < results->capacity; i++) {
            free(results->taxonomies[i]);
        }
        free(results->taxonomies);
        
        omp_destroy_lock(&results->lock);
        free(results);
    }
}

static void add_target_result(TargetResults* results, uint64_t n, const TaxonomyAnalysis* analysis, SearchContext* ctx) {
    omp_set_lock(&results->lock);
    
    // Expand capacity if necessary
    if (results->count >= results->capacity) {
        int old_capacity = results->capacity;
        int new_capacity = results->capacity * 2;
        
        results->n_values = safe_realloc(results->n_values, new_capacity * sizeof(uint64_t), "n values expansion");
        results->taxonomies = safe_realloc(results->taxonomies, new_capacity * sizeof(char*), "taxonomies expansion");
        results->first_mr_pos = safe_realloc(results->first_mr_pos, new_capacity * sizeof(int), "first mr pos expansion");
        results->second_mr_pos = safe_realloc(results->second_mr_pos, new_capacity * sizeof(int), "second mr pos expansion");
        results->max_m_value = safe_realloc(results->max_m_value, new_capacity * sizeof(uint64_t), "max m value expansion");
        results->max_m_position = safe_realloc(results->max_m_position, new_capacity * sizeof(int), "max m position expansion");
        results->max_collatz_value = safe_realloc(results->max_collatz_value, new_capacity * sizeof(uint64_t), "max collatz value expansion");
        results->max_collatz_position = safe_realloc(results->max_collatz_position, new_capacity * sizeof(int), "max collatz position expansion");
        
        // Initialize new taxonomy strings
        for (int i = old_capacity; i < new_capacity; i++) {
            results->taxonomies[i] = safe_malloc(8 * sizeof(char), "taxonomy string");
        }
        
        results->capacity = new_capacity;
    }
    
    results->n_values[results->count] = n;
    results->first_mr_pos[results->count] = analysis->first_mr_pos;
    results->second_mr_pos[results->count] = analysis->second_mr_pos;
    results->max_m_value[results->count] = analysis->max_m_value;
    results->max_m_position[results->count] = analysis->max_m_position;
    results->max_collatz_value[results->count] = analysis->max_collatz_value;
    results->max_collatz_position[results->count] = analysis->max_collatz_position;
    
    // Set taxonomy string
    switch (analysis->taxonomy) {
        case 'A':
            strcpy(results->taxonomies[results->count], "A");
            break;
        case 'B':
            strcpy(results->taxonomies[results->count], "B");
            break;
        case 'C':
            strcpy(results->taxonomies[results->count], "C");
            break;
        default:
            strcpy(results->taxonomies[results->count], "Simple");
            break;
    }
    
    results->count++;
    
    // Only flush to disk if we exceed the threshold AND we have very many results
    // This ensures small and medium datasets stay in memory for detailed display
    if (ctx->progress->found_count >= 50000 && results->count >= RESULTS_FLUSH_INTERVAL) {
        omp_unset_lock(&results->lock);
        flush_results_to_disk(ctx);
    } else {
        omp_unset_lock(&results->lock);
    }
}

static void report_target_found(uint64_t n, uint64_t mr, const TaxonomyAnalysis* analysis, const TargetResults* results, ProgressTracker* tracker) {
    update_last_n_found(tracker, n);
    
    #pragma omp critical(target_report)
    {
        printf(" [*] Found n = %lu: mr = %lu, taxonomy = %c, first_mr_pos = %d, second_mr_pos = %d, max_m = %lu at pos %d (total: %d)\n", 
               n, mr, analysis->taxonomy, analysis->first_mr_pos, analysis->second_mr_pos, 
               analysis->max_m_value, analysis->max_m_position, results->count);
        fflush(stdout);
    }
}

/* ============================================================================
 * COLLATZ SEQUENCE ANALYSIS
 * ============================================================================
 */

static uint64_t find_first_mr_in_sequence(uint64_t n_start, bool* found) {
    uint64_t n = n_start;
    MValues m_values;
    init_m_values(&m_values);
    
    uint64_t first_mr = 0;
    *found = false;
    
    // Generate Collatz sequence and search for repetitions in m values
    for (int step = 0; step < MAX_SEQUENCE_LENGTH && n != 1; step++) {
        uint64_t m = calculate_m(n);
        
        // Check repetition BEFORE adding
        if (is_m_repeated(&m_values, m)) {
            first_mr = m;
            *found = true;
            break;
        }
        
        add_m_value(&m_values, m);
        
        // Apply Collatz transformation with overflow check
        if (!apply_collatz_transform(&n)) {
            *found = false;
            break;
        }
        
        // Safety check for very long sequences
        if (step > 10000 && n > n_start * 1000) {
            first_mr = 0;
            *found = true;
            break;
        }
    }
    
    // If we naturally reach n=1 (trivial cycle), then mr=0
    if (n == 1 && !*found) {
        first_mr = 0;
        *found = true;
    }
    
    destroy_m_values(&m_values);
    return first_mr;
}

static void process_single_number(uint64_t n, SearchContext* ctx, uint64_t* local_found, uint64_t* local_processed) {
    bool mr_found = false;
    uint64_t mr = find_first_mr_in_sequence(n, &mr_found);
    
    if (mr_found && mr == ctx->target_mr) {
        // Get detailed taxonomy analysis for this n
        TaxonomyAnalysis analysis = analyze_sequence_taxonomy(n, ctx->target_mr);
        
        (*local_found)++;
        increment_progress_counters(ctx->progress, true);
        
        // Add to results with detailed analysis
        add_target_result(ctx->results, n, &analysis, ctx);
        
        // Report finding immediately with details (only for small mr values)
        if (ctx->target_mr > 0) {
            report_target_found(n, mr, &analysis, ctx->results, ctx->progress);
        }
    } else {
        increment_progress_counters(ctx->progress, false);
    }
    
    (*local_processed)++;
}

/* ============================================================================
 * PARALLEL SCHEDULING CONFIGURATION
 * ============================================================================ 
 */

static SchedulingConfig configure_parallel_scheduling(uint64_t max_n) {
    SchedulingConfig config;
    config.num_threads = omp_get_max_threads();
    
    // Choose scheduling strategy based on problem size
    if (max_n < 1000000) {
        config.type = SCHEDULE_STATIC;
        config.chunk_size = 0; // Let OpenMP decide
    } else if (max_n < 100000000) {
        config.type = SCHEDULE_GUIDED;
        config.chunk_size = 0; // Not applicable for guided
    } else {
        config.type = SCHEDULE_DYNAMIC;
        // Calculate optimal chunk size
        config.chunk_size = max_n / (config.num_threads * CHUNK_MULTIPLIER);
        if (config.chunk_size < MIN_CHUNK_SIZE) config.chunk_size = MIN_CHUNK_SIZE;
        if (config.chunk_size > MAX_CHUNK_SIZE) config.chunk_size = MAX_CHUNK_SIZE;
    }
    
    return config;
}

static void print_scheduling_info(const SchedulingConfig* config) {
    printf("Using scheduling strategy: ");
    switch(config->type) {
        case SCHEDULE_STATIC:
            printf("static (small range)\n");
            break;
        case SCHEDULE_GUIDED:
            printf("guided (medium range)\n");
            break;
        case SCHEDULE_DYNAMIC:
            printf("dynamic with chunk size %lu (large range)\n", config->chunk_size);
            break;
    }
    printf("Number of threads: %d\n\n", config->num_threads);
}

/* ============================================================================
 * MAIN PARALLEL SEARCH EXECUTION
 * ============================================================================
 */

static void execute_search_with_static_scheduling(SearchContext* ctx, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(static)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_search_with_guided_scheduling(SearchContext* ctx, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(guided)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_search_with_dynamic_scheduling(SearchContext* ctx, uint64_t chunk_size, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_parallel_search(SearchContext* ctx, uint64_t* found_count, uint64_t* processed_count) {
    SchedulingConfig config = configure_parallel_scheduling(ctx->max_n);
    print_scheduling_info(&config);
    
    // Execute search with appropriate scheduling
    switch(config.type) {
        case SCHEDULE_STATIC:
            execute_search_with_static_scheduling(ctx, found_count, processed_count);
            break;
        case SCHEDULE_GUIDED:
            execute_search_with_guided_scheduling(ctx, found_count, processed_count);
            break;
        case SCHEDULE_DYNAMIC:
            execute_search_with_dynamic_scheduling(ctx, config.chunk_size, found_count, processed_count);
            break;
    }
}

/* ============================================================================
 * RESULTS OUTPUT AND REPORTING
 * ============================================================================ 
 */

static int* create_sorting_indices(const TargetResults* results) {
    int* indices = safe_malloc(results->count * sizeof(int), "sorting indices");
    
    for (int i = 0; i < results->count; i++) {
        indices[i] = i;
    }
    
    // Simple bubble sort for n values
    for (int i = 0; i < results->count - 1; i++) {
        for (int j = 0; j < results->count - 1 - i; j++) {
            if (results->n_values[indices[j]] > results->n_values[indices[j + 1]]) {
                int temp = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = temp;
            }
        }
    }
    
    return indices;
}

static void print_detailed_results_from_memory(const SearchContext* ctx) {
    TargetResults* results = ctx->results;
    int* indices = create_sorting_indices(results);
    
    // Count taxonomies
    int taxonomy_counts[4] = {0};
    for (int i = 0; i < results->count; i++) {
        char* tax = results->taxonomies[i];
        if (strcmp(tax, "A") == 0) taxonomy_counts[0]++;
        else if (strcmp(tax, "B") == 0) taxonomy_counts[1]++;
        else if (strcmp(tax, "C") == 0) taxonomy_counts[2]++;
        else taxonomy_counts[3]++;
    }
    
    // Print header for detailed results
    printf(" %8s %8s %9s %10s %10s %11s %14s %16s\n", "n", "Tax", "1st_mr", "2nd_mr", "Max_m", "Max_m_pos", "Max_collatz", "Max_collatz_pos");
    printf(" %8s %8s %9s %10s %10s %11s %14s %16s\n", "--------", "--------", "---------", "----------", "----------", "-----------", "--------------", "----------------");
    
    // Print first 50 results with detailed info, then summarize if more
    int display_limit = (results->count > 50) ? 50 : results->count;
    
    for (int i = 0; i < display_limit; i++) {
        int idx = indices[i];
        printf(" %8lu %8s %9d %10d %10lu %11d %14lu %16d\n", 
               results->n_values[idx], 
               results->taxonomies[idx],
               results->first_mr_pos[idx],
               results->second_mr_pos[idx],
               results->max_m_value[idx],
               results->max_m_position[idx],
               results->max_collatz_value[idx],
               results->max_collatz_position[idx]);
    }
    
    if (results->count > 50) {
        printf(" ... and %d more values\n", results->count - 50);
    }
    
    printf("\nTaxonomy Distribution:\n");
    printf(" Type A: %d (%.1f%%)\n", taxonomy_counts[0], (double)taxonomy_counts[0]/results->count*100);
    printf(" Type B: %d (%.1f%%)\n", taxonomy_counts[1], (double)taxonomy_counts[1]/results->count*100);
    printf(" Type C: %d (%.1f%%)\n", taxonomy_counts[2], (double)taxonomy_counts[2]/results->count*100);
    printf(" Simple: %d (%.1f%%)\n", taxonomy_counts[3], (double)taxonomy_counts[3]/results->count*100);
    
    // Show summary statistics for positions and values
    if (results->count > 0) {
        printf("\nStatistics:\n");
        
        // Calculate averages for non-Simple cases
        int valid_cases = 0;
        double avg_first_pos = 0, avg_second_pos = 0, avg_distance = 0;
        uint64_t min_max_m = UINT64_MAX, max_max_m = 0;
        uint64_t min_max_collatz = UINT64_MAX, max_max_collatz = 0;
        
        for (int i = 0; i < results->count; i++) {
            if (results->first_mr_pos[i] > 0 && results->second_mr_pos[i] > 0) {
                valid_cases++;
                avg_first_pos += results->first_mr_pos[i];
                avg_second_pos += results->second_mr_pos[i];
                avg_distance += (results->second_mr_pos[i] - results->first_mr_pos[i]);
            }
            
            if (results->max_m_value[i] < min_max_m) min_max_m = results->max_m_value[i];
            if (results->max_m_value[i] > max_max_m) max_max_m = results->max_m_value[i];
            
            if (results->max_collatz_value[i] < min_max_collatz) min_max_collatz = results->max_collatz_value[i];
            if (results->max_collatz_value[i] > max_max_collatz) max_max_collatz = results->max_collatz_value[i];
        }
        
        if (valid_cases > 0) {
            avg_first_pos /= valid_cases;
            avg_second_pos /= valid_cases;
            avg_distance /= valid_cases;
            
            printf(" Average first mr position: %.1f\n", avg_first_pos);
            printf(" Average second mr position: %.1f\n", avg_second_pos);
            printf(" Average distance between mr: %.1f\n", avg_distance);
        }
        
        printf(" Smallest maximum m value: %lu\n", min_max_m);
        printf(" Largest maximum m value: %lu\n", max_max_m);
        printf(" Smallest maximum Collatz value: %lu\n", min_max_collatz);
        printf(" Largest maximum Collatz value: %lu\n", max_max_collatz);
    }
    
    // Print comma-separated list of n values
    printf("\nFirst few n values: ");
    int show_count = (results->count > 10) ? 10 : results->count;
    for (int i = 0; i < show_count; i++) {
        int idx = indices[i];
        printf("%lu", results->n_values[idx]);
        if (i < show_count - 1) printf(", ");
    }
    if (results->count > 10) printf(", ...");
    printf("\n");
    
    free(indices);
}

static void print_validation_results(const SearchContext* ctx, double total_time, uint64_t found_count, uint64_t processed_count, const char* filename) {
    printf("\nSEARCH RESULTS:\n");
    printf("\nAll n values that generate mr = %lu:\n\n", ctx->target_mr);
    
    if (found_count > 0) {
        // For small result sets (< 1000), show detailed info from memory
        // For medium result sets (1000-50000), show first 50 with limited stats
        // For very large result sets (> 50000), show simplified info and refer to CSV
        if (found_count < 1000) {
            print_detailed_results_from_memory(ctx);
        } else if (found_count <= 50000) {
            // Show first 50 results even for medium-sized datasets
            printf("Showing first 50 results (total found: %lu):\n\n", found_count);
            if (ctx->results->count > 0) {
                print_detailed_results_from_memory(ctx);
            } else {
                printf("Results have been written incrementally to CSV file.\n");
                printf("Check the output file for complete detailed results grouped by taxonomy.\n\n");
            }
        } else {
            printf("Results have been written incrementally to the CSV file during execution.\n");
            printf("Check the output file for complete detailed results grouped by taxonomy.\n\n");
            printf("Total results found: %lu (dataset too large for detailed display)\n", found_count);
        }
    } else {
        printf("No results found for the specified criteria.\n\n");
    }
    
    printf("\nFinal Summary:\n\n");
    printf(" [*] Target mr value: %lu\n", ctx->target_mr);
    printf(" [*] Search range: 1 to 2^%d = %lu\n", ctx->exponent, processed_count);
    printf(" [*] Total time: %.2f seconds\n", total_time);
    printf(" [*] Speed: %.1f numbers/second\n", (double)processed_count / total_time);
    printf(" [*] Numbers processed: %lu\n", processed_count);
    printf(" [*] Numbers matching criteria: %lu\n", found_count);
    printf(" [*] Percentage matching criteria: %.6f%%\n", (double)found_count / processed_count * 100.0);
}

/* ============================================================================
 * COMMAND LINE ARGUMENT PROCESSING
 * ============================================================================ 
 */

static bool validate_and_parse_arguments(int argc, char* argv[], int* exponent, uint64_t* max_n, uint64_t* target_mr) {
    if (argc != 3) {
        printf("Usage: %s <exponent> <target_mr>\n", argv[0]);
        printf("Examples:\n");
        printf("  %s 25 42     (finds all n < 2^25 where mr = 42)\n", argv[0]);
        printf("  %s 20 0      (finds all n < 2^20 where mr = 0)\n", argv[0]);
        printf("\nRecommended exponents:\n");
        printf("  20 -> 2^20 = 1,048,576 (quick test)\n");
        printf("  25 -> 2^25 = 33,554,432 (default)\n");
        printf("  30 -> 2^30 = 1,073,741,824 (intensive use)\n");
        return false;
    }
    
    *exponent = atoi(argv[1]);
    *target_mr = strtoull(argv[2], NULL, 10);
    
    if (*exponent < 1 || *exponent > 64) {
        printf("Error: Exponent must be between 1 and 64\n");
        printf("Exponent %d is out of valid range\n", *exponent);
        return false;
    }
    
    *max_n = 1UL << *exponent;
    return true;
}

static void print_program_header(int exponent, uint64_t max_n, uint64_t target_mr) {
    printf("\n---------------------------------------------------------------------------------------\n");
    printf("- Searching for n values with specific mr in Collatz sequences                       -\n");
    printf("---------------------------------------------------------------------------------------\n");
    printf("\nUsing %d threads for parallelization\n", omp_get_max_threads());
    printf("Searching for all n < 2^%d that generate mr = %lu\n", exponent, target_mr);
    printf("Range: from 1 to %lu \n\n", max_n - 1);
}

/* ============================================================================
 * MAIN FUNCTION
 * ============================================================================
 */

int main(int argc, char* argv[]) {
    // Parse and validate command line arguments
    int exponent;
    uint64_t max_n, target_mr;
    
    if (!validate_and_parse_arguments(argc, argv, &exponent, &max_n, &target_mr)) {
        return 1;
    }
    
    // Print program header
    print_program_header(exponent, max_n, target_mr);
    
    // Generate output filename BEFORE initializing context
    char filename[256];
    snprintf(filename, sizeof(filename), "values_of_n_lower_than_2_exp_%d_where_mr_%lu.csv", exponent, target_mr);
    
    // Initialize search context
    SearchContext ctx = {
        .max_n = max_n,
        .target_mr = target_mr,
        .exponent = exponent,
        .results = create_target_results(),
        .progress = create_progress_tracker(),
        .start_time = omp_get_wtime(),
        .output_file = NULL
    };
    
    // Initialize incremental file writing
    if (!init_incremental_output(&ctx, filename)) {
        destroy_target_results(ctx.results);
        destroy_progress_tracker(ctx.progress);
        return 1;
    }
    
    // Execute the parallel search
    uint64_t found_count = 0;
    uint64_t processed_count = 0;
    
    execute_parallel_search(&ctx, &found_count, &processed_count);
    
    // For small and medium datasets, write to CSV at the end from memory
    // For very large datasets, results are already written incrementally
    if (found_count < 50000 && ctx.results->count > 0) {
        // Write final results to CSV for small/medium datasets
        write_final_results_to_disk(&ctx);
    }
    
    // Flush any remaining results to disk (for very large datasets)
    flush_results_to_disk(&ctx);
    
    double end_time = omp_get_wtime();
    double total_time = end_time - ctx.start_time;
    
    // Close incremental output
    close_incremental_output(&ctx);
    
    // Print validation results and summary
    print_validation_results(&ctx, total_time, found_count, processed_count, filename);
    
    printf("\n\nResults saved to '%s'\n", filename);
    
    // Cleanup
    destroy_target_results(ctx.results);
    destroy_progress_tracker(ctx.progress);
    
    return 0;
}