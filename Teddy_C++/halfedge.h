/* Based on code written by Kevin Li, Caltech class of 2016 for CS174. */

#pragma once

#include <string>
#include <vector>
#include <math.h>
#include <map>
#include <utility>
#include <algorithm>
#include <cassert>
#include "structs.h"
#include "../Eigen/Dense"


/* Halfedge implementation. */

struct HE // HE for halfedge
{
    // the vertex that this halfedge comes out of
    struct HEV *vertex;
    // the face adjacent to this halfedge (NULL if boundary)
    struct HEF *face;
    // the flip and next halfedge as described in the lecture notes
    struct HE *flip, *next;

};

struct HEF // HEF for halfedge face
{
    // the halfedge associated with this face
    struct HE *he;
    bool visited;
    int type;
};

struct HEV // HEV for halfedge vertex
{
    // the coordinates of the vertex in the mesh
    float x, y, z;
    bool isBoundary;
    // the halfedge going out off this vertex
    struct HE *out;

    int index;
    // store the normal vector for the vertex
    Eigen::Vector3f normal;
};

/* Function prototypes */

static std::pair<int, int> get_edge_key(int x, int y);

static void hash_edge(std::map<std::pair<int, int>, HE*> &edge_hash,
                      std::pair<int, int> edge_key,
                      HE *edge);

static void add_vertex(float x, float y, float z, bool isBoundary, std::vector<HEV*> *hevs);

static void add_face(int idx1, int idx2, int idx3, 
    std::vector<HEV*> *hevs, std::vector<HEF*> *hefs, std::map<std::pair<int, int>, HE*> &edge_hash);

static void delete_HEV(std::vector<HEV*> *hevs);
static void delete_HEF(std::vector<HEF*> *hefs);



/* Function implementations */

/* Order the indices of the key so mapping can be unique. */
static std::pair<int, int> get_edge_key(int x, int y)
{
    assert(x != y);
    return std::pair<int, int>(std::min(x, y), std::max(x, y));
}

/* Map from two vertices (idx1, idx2) to its associated halfedges. */
static void hash_edge(std::map<std::pair<int, int>, HE*> &edge_hash,
                     std::pair<int, int> edge_key,
                     HE *edge)
{   
    // Find the halfedge that already lives at the given edge.
    if(edge_hash.count(edge_key) > 0)
    {
        HE *flip = edge_hash[edge_key];
        flip->flip = edge;
        edge->flip = flip;
    }
    // If this edge doesn't have an associated halfedge yet, add the current one.
    else {
        edge_hash[edge_key] = edge;
    }

}


static void add_vertex(float x, float y, float z, bool isBoundary, std::vector<HEV*> *hevs) {

    HEV *hev = new HEV;
    hev->x = x;
    hev->y = y;
    hev->z = z;
    hev->isBoundary = isBoundary;
    hev->out = NULL;
    hev->index = hevs->size();
    hevs->push_back(hev);
}

static void add_face(int idx1, int idx2, int idx3, 
    std::vector<HEV*> *hevs, std::vector<HEF*> *hefs, 
    std::map<std::pair<int, int>, HE*> &edge_hash) {

    HE *e1 = new HE;
    HE *e2 = new HE;
    HE *e3 = new HE;

    HEF *f = new HEF;
    f->visited = false;
    f->he = e1;

    e1->face = f;
    e2->face = f;
    e3->face = f;

    e1->flip = NULL;
    e2->flip = NULL;
    e3->flip = NULL;

    e1->next = e2;
    e2->next = e3;
    e3->next = e1;

    e1->vertex = hevs->at(idx1);
    e2->vertex = hevs->at(idx2);
    e3->vertex = hevs->at(idx3);

    hevs->at(idx1)->out = e1;
    hevs->at(idx2)->out = e2;
    hevs->at(idx3)->out = e3;

    hash_edge(edge_hash, get_edge_key(idx1, idx2), e1);
    hash_edge(edge_hash, get_edge_key(idx2, idx3), e2);
    hash_edge(edge_hash, get_edge_key(idx3, idx1), e3);

    hefs->push_back(f);
}

static void delete_HEV(std::vector<HEV*> *hevs)
{
    int hev_size = hevs->size();
    

    for(int i = 1; i < hev_size; ++i) {
        delete hevs->at(i);
    }

    delete hevs;
    
}

static void delete_HEF(std::vector<HEF*> *hefs) {

    int num_hefs = hefs->size();

    for(int i = 0; i < num_hefs; ++i)
    {
        delete hefs->at(i)->he->next->next;
        delete hefs->at(i)->he->next;
        delete hefs->at(i)->he;
        delete hefs->at(i);
    }
    delete hefs;
}

