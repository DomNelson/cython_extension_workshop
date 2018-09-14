#ifndef SIGNAL
#define SIGNAL
#include <gsl/gsl_rng.h>

typedef struct node_t_t {
    int ID;
    double weight;
    int genotype;

    struct node_t_t *mother;
    struct node_t_t *father;
    int climbed_to_mother;
    int climbed_to_father;

    int *active_samples; // Store length of array here as well as ped?
} node_t;

typedef struct {
    int idx;
    node_t *node;
    char status; // A - active, C - coalesced, F - founder
} lineage_t;

// Not used/needed?
typedef struct sample_list_t_t {
    node_t *node;
    struct sample_list_t_t *next;
} sample_list_t;

typedef struct {
    uint32_t num_nodes;
    uint32_t num_samples;
    uint32_t num_active_lineages;
    double sim_homs;
    node_t *node_array;
    node_t **samples;
    lineage_t *active_lineages;
    gsl_rng *rng;
} ped_t;

void multiply_by_10_in_C(double arr[], unsigned int n);

ped_t *ped_alloc(void);
int ped_nodes_alloc(ped_t *ped, uint32_t num_nodes, int num_samples);
int ped_samples_alloc(ped_t *ped, uint32_t num_samples);
int ped_alloc_rng(ped_t *ped);
int free_ped(ped_t *ped);

int ped_load(ped_t *ped, int *inds, int *fathers, int *mothers, int num_inds);
int ped_load_samples_from_idx(ped_t *ped, int *samples_idx, int *genotypes, int num_samples);
int ped_init_sample_weights(ped_t *ped);

int ped_print_nodes(ped_t *ped);
int ped_print_samples(ped_t *ped);
void print_node(node_t node);

int ped_update_ancestor_weights(ped_t *ped, node_t *node, int sample_idx, double delta);
int ped_update_ancestor_weights_from_idx(ped_t *ped, int node_idx, double delta);
int ped_set_all_weights(ped_t *ped, double val);
double ped_get_node_weight_from_idx(ped_t *ped, int node_idx);
double node_get_parent_weight(node_t *node, int gen);
int update_parent_carrier(ped_t *ped, node_t *node, int sample_idx);
int update_parent_not_carrier(ped_t *ped, node_t *node, int sample_idx);
int update_parent_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx);
int update_parent_not_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx);
int node_get_max_coalescences(ped_t *ped, node_t *node);

int ped_climb_step(ped_t *ped);
int ped_lineage_coalesce(ped_t *ped, lineage_t *lineage);
int ped_lineage_climb(ped_t *ped, lineage_t *lineage);
int ped_lineage_set_next_parent(ped_t *ped, lineage_t *lineage, char parent);
int ped_lineage_update_genotype_founder(ped_t *ped, lineage_t *lineage);
#endif
