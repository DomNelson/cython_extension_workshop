#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_minmax.h>

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
    ped->num_active_lineages = 0;
    ped->node_array = NULL;
    ped->samples = NULL;
    ped->active_lineages = NULL;

    return ped;
}

int node_init(node_t *node, int num_samples) {
    node->ID = -1;
    node->weight = 0;
    node->genotype = 0;
    node->climbed_to_mother = 0;
    node->climbed_to_father = 0;
    node->father = NULL;
    node->mother = NULL;
    node->active_samples = calloc(num_samples, sizeof(int));

    return 0;
}

int free_ped(ped_t *ped) {
    int i;
    int count = 0;
    node_t *node;

    if (ped->samples != NULL) {
        printf("Freeing ped->samples\n");
        free(ped->samples);
    }
    if (ped->active_lineages != NULL) {
        printf("Freeing ped->active_lineages\n");
        free(ped->active_lineages);
    }
    for (i = 0; i < ped->num_nodes; i++) {
        node = &ped->node_array[i];

        if (node->active_samples != NULL) {
            count++;
            free(node->active_samples);
        }
    }
    if (count > 0) {
        printf("Freed %d node->active_samples\n", count);
    }
    if (ped->node_array != NULL) {
        printf("Freeing ped->node_array\n");
        free(ped->node_array);
    }
    gsl_rng_free(ped->rng);
    free(ped);

    return 0;
}

int ped_nodes_alloc(ped_t *ped, uint32_t num_nodes, int num_samples) {
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
        node_init(n, num_samples);
    }
out:
    return ret;
}

int ped_samples_alloc(ped_t *ped, uint32_t num_samples) {
    int ret = 0;
    int i;

    ped->num_samples = num_samples;
    ped->num_active_lineages = num_samples;

    ped->samples = calloc(num_samples, sizeof(node_t*));
    if (ped->samples == NULL){
        ret = 1;
        goto out;
    }
    ped->active_lineages = calloc(num_samples, sizeof(lineage_t));
    if (ped->active_lineages == NULL){
        ret = 1;
        goto out;
    }
    printf("Allocated active_lineages\n");
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
        n->ID = ID;
        n->father = father;
        n->mother = mother;
    }

    return ret;
}

int ped_load_samples_from_idx(ped_t *ped, int *samples_idx, int *genotypes, int num_samples) {
    int ret = 0;
    int i, s_idx;
    node_t *n = NULL;
    lineage_t *l = NULL;

    assert(num_samples == ped->num_samples);
    printf("Loading %d samples\n", num_samples);

    for (i = 0; i < num_samples; i++) {
        s_idx = samples_idx[i];
        n = &ped->node_array[s_idx];
        ped->samples[i] = n;
        printf("%d\n", i);

        l = &ped->active_lineages[i];
        l->idx = s_idx;
        l->node = n;
        l->status = 'A'; // Initial state is 'active'

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

int ped_print_lineages(ped_t *ped) {
    int ret = 0;
    int i, j;
    lineage_t *l;
    node_t *n;

    printf("Active\n");
    for (i = 0; i < ped->num_samples; i++)  {
        if (i == ped->num_active_lineages) {
            printf("Inactive\n");
        }
        l = &ped->active_lineages[i];
        n = l->node;

        printf("Lineage idx: %d, node ID: %d, status: %c\n",
                l->idx, n->ID, l->status);
    }
    return ret;
}

int ped_print_nodes(ped_t *ped) {
    int ret = 0;
    int i, j;
    node_t *n;

    for (i = 0; i < ped->num_nodes; i++)  {
        n = &ped->node_array[i];
        printf("Node ID: %d, weight: %f, genotype: %d",
                n->ID, n->weight, n->genotype);

        if (n->active_samples != NULL) {
            printf(", active_samples: [ ");
            for (j = 0; j < ped->num_samples; j++) {
                printf("%d ", n->active_samples[j]);
            }
            printf("]");
        } else {
            printf(", active_samples: not allocated");
        }

        printf(", max_coal %d", node_get_max_coalescences(ped, n));

        if (n->father != NULL) {
            printf(", father: %d", n->father->ID);
        }
        if (n->mother != NULL) {
            printf(", mother: %d\n", n->mother->ID);
        } else {
            printf(", founder\n");
        }
    }
    ped_print_lineages(ped);
    return ret;
}

node_t* ped_get_node(ped_t *ped, int node_idx) {
    node_t *node;

    node = &ped->node_array[node_idx];

    return node;
}

int ped_update_ancestor_weights(ped_t *ped, node_t *node, int sample_idx, double delta) {
    // TODO: This will eventually need to be updated - currently only tracks a
    // single index per lineage, but in fact we need to track the index of this
    // lineage as well as all those which have coalesced with it
    // --> alternatively, update active_samples of all nodes to 0 when a
    // lineage coalesces, and decrement minimum number of possible coalescences

    int ret = 0;
    node_t *mother, *father;

    /* printf("Updating ancestor %d, with index %d\n", */
    /*         node->ID, sample_idx); */

    assert(delta != 0);
    assert(sample_idx >= 0);
    assert(node->active_samples != NULL);
    node->weight = node->weight + delta;

    if (delta < 0) {
        assert(node->active_samples[sample_idx] == 1);
        node->active_samples[sample_idx] = 0;
    } else if (delta > 0) {
        node->active_samples[sample_idx] = 1;
    }

    if (node->father != NULL) {
        ret = ped_update_ancestor_weights(ped, node->father, sample_idx, delta / 2);
    }
    if (node->mother != NULL) {
        ret = ped_update_ancestor_weights(ped, node->mother, sample_idx, delta / 2);
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

int ped_init_sample_weights(ped_t *ped) {
    int ret = 0;
    int i;
    lineage_t *lineage;

    for (i = 0; i < ped->num_active_lineages; i++) {
        lineage = &ped->active_lineages[i];
        lineage->idx= i;
        ret = ped_update_ancestor_weights(ped, lineage->node, lineage->idx, 1);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}

double node_get_parent_weight(node_t *node, int gen) {
    int ret = 0;
    double weight;
    node_t *mother, *father;

    assert(node != NULL);

    // Since this is a parent, the weight includes signal from the offspring
    // being climbed, which we must subtract. Subsequent generations have this
    // subtraction adjusted automatically by the generation coefficient.
    weight = node->weight - 0.5;

    if (node->father != NULL) {
        weight += pow(2, -gen) * node_get_parent_weight(node->father, gen + 1);
    }

    if (node->mother != NULL) {
        weight += pow(2, -gen) * node_get_parent_weight(node->mother, gen + 1);
    }

    return weight;
}

double ped_get_node_weight_from_idx(ped_t *ped, int node_idx) {
    double weight = 0;
    node_t *node;

    node = ped_get_node(ped, node_idx);
    weight = node_get_parent_weight(node, 0);

    return weight;
}

int update_parent_carrier(ped_t *ped, node_t *node, int sample_idx) {
    int ret = 0;

    ret = ped_update_ancestor_weights(ped, node, sample_idx, 0.5);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int update_parent_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx) {
    int ret = 0;
    node_t *node;

    node = ped_get_node(ped, node_idx);
    ret = ped_update_ancestor_weights(ped, node, sample_idx, 0.5);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int update_parent_not_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx) {
    int ret = 0;
    node_t *node;

    node = ped_get_node(ped, node_idx);
    ret = ped_update_ancestor_weights(ped, node, sample_idx, -0.5);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int update_parent_not_carrier(ped_t *ped, node_t *node, int sample_idx) {
    int ret = 0;

    ret = ped_update_ancestor_weights(ped, node, sample_idx, -0.5);
    if (ret != 0) {
        goto out;
    }
out:
    return ret;
}

int node_get_max_coalescences(ped_t *ped, node_t *node) {
    int i;
    int max_coal;

    assert(node != NULL);
    max_coal = 0;
    /* printf("Summing %d possible coalescences\n", ped->num_samples); */
    for (i = 0; i < ped->num_samples; i++) {
        /* printf("%d\n", i); */
        max_coal = max_coal + node->active_samples[i];
    }

    if (node->father != NULL) {
        max_coal = GSL_MAX_INT(
                max_coal, node_get_max_coalescences(ped, node->father));
    }

    if (node->mother != NULL) {
        max_coal = GSL_MAX_INT(
                max_coal, node_get_max_coalescences(ped, node->mother));
    }

    return max_coal;
}

int ped_lineage_update_genotype(ped_t *ped, lineage_t *lineage) {
    // TODO: Pass IS factor as arg to update?
    int ret = 0;
    node_t *node;
    int mother_num_coal, father_num_coal;
    double loglik = 1; // use ped->loglik?

    node = lineage->node;
    if (node->mother == NULL && node->father == NULL) {
        ret = ped_lineage_update_genotype_founder(ped, lineage);
        goto out;
    }

    if (node->genotype == 0) {
        node->genotype = 1;
    } else if (node->genotype == 1) {
        // Need to check if we can create a homozygote
        printf("Checking if we can make a homozygote in %d\n",
                node->ID);
        mother_num_coal = father_num_coal = 0;
        if (node->mother != NULL) {
            mother_num_coal = node_get_max_coalescences(ped, node->mother);
        }
        if (node->father != NULL) {
            father_num_coal = node_get_max_coalescences(ped, node->father);
        }

        if (mother_num_coal < ped->num_samples ||
                father_num_coal < ped->num_samples) {
            // Can't create a homozygote - coalesce and we're done
            loglik = loglik * 0.5;
            ped_lineage_coalesce(ped, lineage);
        } else {
            if (gsl_rng_uniform(ped->rng) < ped->sim_homs) {
                // Homozygote sampled
                loglik = loglik * (0.5 / ped->sim_homs);
                node->genotype = 2;
                // Now climb this lineage again to update possible coalescence
                // points for other lineages
                if (node->climbed_to_mother == 1) {
                    ped_lineage_set_next_parent(ped, lineage, 'f');
                } else if (node->climbed_to_father == 1) {
                    ped_lineage_set_next_parent(ped, lineage, 'm');
                } else {
                    printf("Potential homozygote - other allele hasn't climbed yet\n");
                    // The other lineage hasn't climbed yet. Since one
                    // lineage has to go each way, we can choose a parent
                    // uniformly regardless of weights
                    if (gsl_rng_uniform(ped->rng) < 0.5) {
                        ped_lineage_set_next_parent(ped, lineage, 'm');
                    } else {
                        ped_lineage_set_next_parent(ped, lineage, 'm');
                    }
                }
            } else {
                // No homozygote sampled - coalesce
                loglik = loglik * (0.5 / (1 - ped->sim_homs));
                ped_lineage_coalesce(ped, lineage);
            }
        }
    } else if (node->genotype == 2) {
        ped_lineage_coalesce(ped, lineage);
    }
out:
    return ret;
}

int ped_lineage_update_genotype_founder(ped_t *ped, lineage_t *lineage) {
    int ret = 0;
    node_t *node;
    double loglik = 1;

    node = lineage->node;
    assert(node->mother == NULL && node->father == NULL);

    if (node->genotype == 0) {
        node->genotype = 1;
    } else if (node->genotype == 1) {
        if (gsl_rng_uniform(ped->rng) < ped->sim_homs) {
            // Homozygote sampled
            loglik = loglik * (0.5 / ped->sim_homs);
            node->genotype = 2;
        } else {
            // No homozygote sampled - coalesce
            loglik = loglik * (0.5 / (1 - ped->sim_homs));
            ped_lineage_coalesce(ped, lineage);
        }
    } else if (node->genotype == 2) {
        ped_lineage_coalesce(ped, lineage);
    }

    //TODO: Think about whether we should call reached_founder
    // here or not...

    return ret;
}

int ped_lineage_reached_founder(ped_t *ped, lineage_t *lineage) {
    int ret = 0;
    lineage_t tmp;
    lineage_t *last;

    // Swap current lineage with last active lineage
    last = ped->active_lineages + ped->num_active_lineages - 1;
    tmp = *last;
    *last = *lineage;
    *lineage = tmp;

    // Reduce number of active lineages by 1 and update status
    assert(lineage->node->mother == NULL && lineage->node->father == NULL);
    printf("---Lineage %d reached founder %d\n",
            lineage->idx, lineage->node->ID);
    ped->num_active_lineages--;
    lineage->status = 'F';

    return ret;
}

int ped_lineage_coalesce(ped_t *ped, lineage_t *lineage) {
    int ret = 0;
    lineage_t tmp;
    lineage_t *last;

    // Swap current lineage with last active lineage
    last = ped->active_lineages + ped->num_active_lineages - 1;
    tmp = *last;
    *last = *lineage;
    *lineage = tmp;

    // Reduce number of active lineages by 1 and update status
    printf("Coalescing lineage %d in node %d\n",
            last->idx, last->node->ID);
    ped->num_active_lineages--;
    last->status = 'C';

    return ret;
}

int ped_lineage_set_next_parent(ped_t *ped, lineage_t *lineage, char parent) {
    int ret = 0;
    node_t *node;

    node = lineage->node;

    if (parent == 'f') {
        assert(node->climbed_to_father == 0);
        update_parent_not_carrier(ped, node->mother, lineage->idx);
        update_parent_carrier(ped, node->father, lineage->idx);
        node->climbed_to_father = 1;
        lineage->node = node->father;
    } else if (parent == 'm') {
        assert(node->climbed_to_mother == 0);
        update_parent_not_carrier(ped, node->father, lineage->idx);
        update_parent_carrier(ped, node->mother, lineage->idx);
        node->climbed_to_mother = 1;
        lineage->node = node->mother;
    } else {
        printf("Error - incorrectly specified parent type: %c\n",
                parent);
        assert(1 == 0);
    }
    ped_lineage_update_genotype(ped, lineage);

    return ret;
}

int ped_lineage_climb(ped_t *ped, lineage_t *lineage) {
    int ret = 0;
    node_t *node;
    int sample_idx;
    int mother_num_coal, father_num_coal;
    double mother_weight, father_weight;
    double x;

    node = lineage->node;
    assert(node != NULL);

    if (node->mother == NULL && node->father == NULL) {
        ped_lineage_reached_founder(ped, lineage);
        goto out;
    }

    mother_weight = father_weight = 0;
    mother_num_coal = father_num_coal = 0;
    if (node->mother != NULL) {
        mother_weight = node_get_parent_weight(node->mother, 0);
        mother_num_coal = node_get_max_coalescences(ped, node->mother);
    }
    if (node->father != NULL) {
        father_weight = node_get_parent_weight(node->father, 0);
        father_num_coal = node_get_max_coalescences(ped, node->father);
    }

    // TODO: Need better handling of when mother/father is NULL
    printf("%d choosing from %d: %f or %d %f\n",
            node->ID,
            node->mother->ID, mother_weight,
            node->father->ID, father_weight);

    assert(mother_weight + father_weight > 0);

    // These options are based purely on laws of inheritance - no IS needed
    if (node->climbed_to_mother == 1) {
        assert(father_num_coal == ped->num_samples);
        ped_lineage_set_next_parent(ped, lineage, 'f');
        goto out;
    }
    if (node->climbed_to_father == 1) {
        assert(mother_num_coal == ped->num_samples);
        ped_lineage_set_next_parent(ped, lineage, 'm');
        goto out;
    }

    // TODO: Add IS factors
    if (mother_num_coal < ped->num_samples) {
        assert(father_num_coal == ped->num_samples);
        ped_lineage_set_next_parent(ped, lineage, 'f');
        goto out;
    }
    if (father_num_coal < ped->num_samples) {
        assert(mother_num_coal == ped->num_samples);
        ped_lineage_set_next_parent(ped, lineage, 'm');
        goto out;
    }

    x = gsl_rng_uniform(ped->rng);
    if (x < mother_weight / (mother_weight + father_weight)) {
        ped_lineage_set_next_parent(ped, lineage, 'm');
    } else {
        ped_lineage_set_next_parent(ped, lineage, 'f');
        goto out;
    }
out:
    return ret;
}

int ped_climb_step(ped_t *ped) {
    int ret = 0;
    int i, j, n;
    lineage_t tmp;
    lineage_t *lineage;

    // Swap random node with last node and choose from remaining.
    // Repeat until done

    // Use local variable for loop since num_active_lineages changes
    // if a lineage coalesces
    n = ped->num_active_lineages;
    for (i = n - 1; i >= 0; i--) {
        j = gsl_rng_uniform_int(ped->rng, i + 1);
        tmp = ped->active_lineages[j];
        ped->active_lineages[j] = ped->active_lineages[i];
        ped->active_lineages[i] = tmp;

        printf("climbing lineage %d array index %d\n",
                ped->active_lineages[i].idx, i);
        assert(ped->active_lineages[i].status == 'A');
        ret = ped_lineage_climb(ped, &ped->active_lineages[i]);
        if (ret != 0) {
            goto out;
        }
    }
out:
    return ret;
}
