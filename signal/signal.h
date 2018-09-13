#ifndef SIGNAL
#define SIGNAL
#include <gsl/gsl_rng.h>

typedef struct node_t_t {
    uint32_t ID;
    double weight;
    uint32_t genotype;
    struct node_t_t *mother;
    struct node_t_t *father;
} node_t;

// Not used/needed?
typedef struct sample_list_t_t {
    node_t *node;
    struct sample_list_t_t *next;
} sample_list_t;

typedef struct {
    uint32_t num_nodes;
    uint32_t num_samples;
    uint32_t num_active_nodes;
    node_t *node_array;
    node_t **samples;
    node_t **active_nodes;
    gsl_rng *rng;
} ped_t;

void multiply_by_10_in_C(double arr[], unsigned int n);

ped_t *ped_alloc(void);
int ped_nodes_alloc(ped_t *ped, uint32_t num_nodes);
int ped_samples_alloc(ped_t *ped, uint32_t num_samples);
int ped_alloc_rng(ped_t *ped);
int free_ped(ped_t *ped);

int ped_load(ped_t *ped, int *inds, int *fathers, int *mothers, int num_inds);
int ped_load_samples_from_idx(ped_t *ped, int *samples_idx, int *genotypes, int num_samples);

int ped_print_nodes(ped_t *ped);
int ped_print_samples(ped_t *ped);
void print_node(node_t node);

int ped_update_ancestor_weights(ped_t *ped, node_t *node, double delta);
int ped_update_ancestor_weights_from_idx(ped_t *ped, int node_idx, double delta);
int ped_set_all_weights(ped_t *ped, double val);
double ped_get_node_weight_from_idx(ped_t *ped, int node_idx);
double node_get_weight(node_t *node, int gen);
int update_parent_carrier_from_idx(ped_t *ped, int node_idx);
int update_parent_not_carrier_from_idx(ped_t *ped, int node_idx);

int ped_climb_step(ped_t *ped);
#endif
