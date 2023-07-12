#include<bits/stdc++.h>
#define pb push_back
using namespace std;

struct HalfEdge;

struct Vertex{
    Vertex(){};
    // Vertex(pair<double, double> p, int k) {
    //     x = p.first;
    //     y = p.second;
    //     key = k;
    // }
    Vertex(const Vertex &v1);
    // Vertex(const Vertex &v1){
    //     for(int i = 0; i<v1.incidentEdges.size(); i++){
    //         HalfEdge *temp = new HalfEdge;
    //         *temp = *v1.incidentEdges[i];
    //         incidentEdges.pb(temp);
    //     }
    //     x = v1.x;
    //     y = v1.y;
    //     // *incidentEdges = v1.incidentEdges;
    //     key = v1.key;
    // }
    double x, y;
    // struct Vertex* pr = NULL;
    int key;
    // struct Vertex* ne = NULL;

    vector<HalfEdge*> incidentEdges;
    Vertex* next;
};

struct Face{
    Face(){};
    Face(const Face &f1);
    // Face(const Face &f1){
    //     edge = new HalfEdge;
    //     *edge = *f1.edge;
    // }
    HalfEdge* edge;
};

struct HalfEdge{
    HalfEdge(){};
    // HalfEdge(Vertex* v) {
    //     origin = v;
    //     origin_vkey = vkey;
    // }
    HalfEdge(const HalfEdge &hf);
    // HalfEdge(const HalfEdge &hf){
    //     origin = new Vertex;
    //     *origin = *(hf.origin);
    //     twin = new HalfEdge;
    //     *twin = *(hf.twin);
    //     next = new HalfEdge;
    //     *next = *hf.next;
    //     prev = new HalfEdge;
    //     *prev = *hf.prev;
    //     face = new Face;
    //     *face = *hf.face;
    //     // origin_vkey = new int;
    //     origin_vkey = hf.origin_vkey;
    // }
    Vertex *origin;
    HalfEdge* twin = NULL;
    HalfEdge* next = NULL;
    HalfEdge* prev = NULL;
    Face* face = NULL;
    int origin_vkey;
};

Vertex::Vertex(const Vertex &v1){
        for(int i = 0; i<v1.incidentEdges.size(); i++){
            HalfEdge *temp = new HalfEdge;
            *temp = *v1.incidentEdges[i];
            incidentEdges.pb(temp);
        }
        double *temp = new double;
        *temp = v1.x; 
        x = *temp;

        temp = new double;
        *temp = v1.y; 
        y = *temp;

        int *temp1 = new int;
        *temp1 = v1.key;
        // *incidentEdges = v1.incidentEdges;
        key = *temp1;
    }

HalfEdge::HalfEdge(const HalfEdge &hf){
        origin = new Vertex;
        *origin = *(hf.origin);
        twin = new HalfEdge;
        *twin = *(hf.twin);
        next = new HalfEdge;
        *next = *hf.next;
        prev = new HalfEdge;
        *prev = *hf.prev;
        face = new Face;
        *face = *hf.face;
        // origin_vkey = new int;
        origin_vkey = hf.origin_vkey;
    }

Face::Face(const Face &f1){
        edge = new HalfEdge;
        *edge = *f1.edge;
    }

struct DCEL{
    DCEL(){};
    DCEL(const DCEL &d){
        for(int i = 0; i < d.vertices.size(); i++){
            Vertex *temp = new Vertex;
            *temp = *d.vertices[i];
            vertices.pb(temp);
        }

        for(int i = 0; i < d.halfEdges.size(); i++){
            HalfEdge *temp = new HalfEdge;
            *temp = *d.halfEdges[i];
            halfEdges.pb(temp);
        }

        for(int i = 0; i < d.faces.size(); i++){
            Face *temp = new Face;
            *temp = *d.faces[i];
            faces.pb(temp);
        }
    }
    vector<Vertex*> vertices;
    vector<HalfEdge*> halfEdges;
    vector<Face*> faces;
    vector<HalfEdge*> diags;

    // void remove_vertex(Vertex v1){
    //     HalfEdge *temp, *temptwin;
    //     for(auto vert: vertices){
    //         if(vert->x == v1.x && vert->y == v1.y){
    //             temp = vert->incidentEdges[1];
    //             temptwin = vert->incidentEdges[0];
    //             break;
    //         }
    //     }
        
    // }
};
