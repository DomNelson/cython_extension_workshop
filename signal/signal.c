#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "signal.h"

int ped_alloc_rng(ped_t *ped) {
    int ret = 0;
    struct timespec ts;
    gsl_rng * r;

    r = gsl_rng_alloc(gsl_rng_taus);
    if (r == NULL) {
        ret = 1;
        goto out;
    }
    clock_gettime(CLOCK_MONOTONIC, &ts);
    gsl_rng_set(r, (time_t)ts.tv_nsec);
    printf("Random seed: %ld\n", (time_t)ts.tv_nsec);

    ped->rng = r;
out:
    return ret;
}

void multiply_by_10_in_C(double arr[], unsigned int n)
{
    unsigned int i;
    for (i = 0; i < n; i++) {
        arr[i] *= 10;
    }
}

void print_node(node_t node) {
    printf("%d\n", node.ID);
    printf("%f\n", node.weight);
}

ped_t * ped_alloc(void) {
    int ret = 0;
    static ped_t *ped;

    ped = calloc(1, sizeof(ped_t));
    assert(ped != NULL);
    ret = ped_alloc_rng(ped);
    assert(ret == 0);


    // Maybe all these should be calloc'd and freed as well?
    ped->num_nodes = 0;
    ped->num_samples = 0;
    ped->num_active_nodes = 0;
    ped->node_array = NULL;
    ped->samples = NULL;
    ped->active_nodes = NULL;

    return ped;
}

int set_node(node_t *node, uint32_t ID, double weight, uint32_t genotype, node_t *mother, node_t *father) {
    node->ID = ID;
    node->weight = weight;
    node->father = father;
    node->mother = mother;

    return 0;
}

int free_ped(ped_t *ped) {
    if (ped->node_array != NULL) {
        printf("Freeing ped->node_array\n");
        free(ped->node_array);
    }
    if (ped->samples != NULL) {
        printf("Freeing ped->samples\n");
        free(ped->samples);
    }
    if (ped->active_nodes != NULL) {
        printf("Freeing ped->active_nodes\n");
        free(ped->active_nodes);
    }
    gsl_rng_free(ped->rng);
    free(ped);

    return 0;
}

int ped_nodes_alloc(ped_t *ped, uint32_t num_nodes) {
    int ret = 0;
    int i;
    node_t *a, *n;

    ped->num_nodes = num_nodes;
    ped->node_array = calloc(num_nodes, sizeof(node_t));
    if (ped->node_array == NULL){
        ret = 1;
        goto out;
    }

    i = 0;
    for (a = ped->node_array; i < num_nodes; i++) {
        n = a + i;
        set_node(n, 0, 0, 0, NULL, NULL);
    }
out:
    return ret;
}

int ped_samples_alloc(ped_t *ped, uint32_t num_samples) {
    int ret = 0;
    int i;

    ped->num_samples = num_samples;
    ped->num_active_nodes = num_samples;

    ped->samples = calloc(num_samples, sizeof(node_t*));
    if (ped->samples == NULL){
        ret = 1;
        goto out;
    }
    ped->active_nodes = calloc(num_samples, sizeof(node_t*));
    if (ped->active_nodes == NULL){
        ret = 1;
        goto out;
    }
    printf("Allocated active_nodes\n");
out:
    return ret;
}

int ped_load(ped_t *ped, int *inds, int *fathers, int *mothers, int num_inds) {
    int ret = 0;
    int i;
    node_t *n;
    int ID;
    node_t *father, *mother;

    for (i = 0; i < num_inds; i++) {
        n = &ped->node_array[i];
        ID = inds[i];
        father = NULL;
        mother = NULL;

        if (fathers[i] != -1) {
            father = &ped->node_array[fathers[i]];
        } else {
        }
        if (mothers[i] != -1) {
            mother = &ped->node_array[mothers[i]];
        }

        set_node(n, ID, 0, 0, father, mother);
    }

    return ret;
}

int ped_load_samples_from_idx(ped_t *ped, int *samples_idx, int *genotypes, int num_samples) {
    int ret = 0;
    int i, s_idx;
    node_t *n;

    assert(num_samples == ped->num_samples);
    printf("Loading %d samples\n", num_samples);

    for (i = 0; i < num_samples; i++) {
        s_idx = samples_idx[i];
        n = &ped->node_array[s_idx];
        ped->samples[i] = n;
        ped->active_nodes[i] = n;

        n->genotype = genotypes[i];
    }
    printf("Done loading samples\n");
    return ret;
}

int ped_print_samples(ped_t *ped) {
    int ret = 0;
    int i;
    node_t n;

    for (i = 0; i < ped->num_samples; i++)  {
        n = *ped->samples[i];
        printf("Sample ID: %d, weight: %f\n", n.ID, n.weight);
    }
    return ret;
}

int ped_print_nodes(ped_t *ped) {
    int ret = 0;
    int i;
    node_t n;

    for (i = 0; i < ped->num_nodes; i++)  {
        n = ped->node_array[i];
        printf("Node ID: %d, weight: %f, genotype: %d", n.ID, n.weight, n.genotype);

        if (n.father != NULL) {
            printf(", father: %d", n.father->ID);
        }
        if (n.mother != NULL) {
            printf(", mother: %d\n", n.mother->ID);
        } else {
            printf("\n");
        }
    }
    return ret;
}

node_t* ped_get_node(ped_t *ped, int node_idx) {
    node_t *node;

    node = &ped->node_array[node_idx];

    return node;
}

int ped_update_ancestor_weights_from_idx(ped_t *ped, int node_idx, double delta) {
    int ret = 0;
    node_t *node;

    node = ped_get_node(ped, node_idx);
    ret = ped_update_ancestor_weights(ped, node, delta);

    return ret;
}

int ped_update_ancestor_weights(ped_t *ped, node_t *node, double delta) {
    int ret = 0;
    node_t *mother, *father;

    /* printf("ID: %d\n", node->ID); */
    /* printf("Old weight: %f\n", node->weight); */
    node->weight = node->weight + delta;
    /* printf("New weight: %f\n\n", node->weight); */

    if (node->father != NULL) {
        ret = ped_update_ancestor_weights(ped, node->father, delta / 2);
    }

    if (node->mother != NULL) {
        ret = ped_update_ancestor_weights(ped, node->mother, delta / 2);
    }

    return ret;
}

int ped_set_all_weights(ped_t *ped, double val) {
    int ret = 0;
    int i;

    for (i = 0; i < ped->num_nodes; i++) {
        ped->node_array[i].weight = val;
    }

    return ret;
}

double ped_get_node_weight_from_idx(ped_t *ped, int node_idx) {
    double weight = 0;
    node_t *node;

    node = ped_get_node(ped, node_idx);
    weight = node_get_weight(node, 0);

    return weight;
}

double node_get_weight(node_t *node, int gen) {
    int ret = 0;
    double weight;
    node_t *mother, *father;

    /* printf("ID: %d\n", node->ID); */
    /* printf("Adding weight: %f\n", pow(2, -gen) * node->weight); */
    weight = node->weight;

    if (node->father != NULL) {
        weight += pow(2, -gen) * node_get_weight(node->father, gen + 1);
    }

    if (node->mother != NULL) {
        weight += pow(2, -gen) * node_get_weight(node->mother, gen + 1);
    }

    return weight;
}
