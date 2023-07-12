#include <bits/stdc++.h>
using namespace std;
#include "dcel_hvg.hpp"
#define pb push_back


bool vertexEqual(Vertex v1, Vertex v2){
    if(v1.x == v2.x && v1.y == v2.y && v1.key == v2.key) // check for v1.key == v2.key
        return true;

    return false;
}

double angle_calc(Vertex a, Vertex b, Vertex c)
{   

    // remove comments
	// Stores coefficient of X
	// direction of vector A[1]A[0]
	// int X1 = (A[1][0] - A[0][0]);
    double x1 = b.x - a.x;
	// Stores coefficient of Y
	// direction of vector A[1]A[0]
	// int Y1 = (A[1][1] - A[0][1]);
    double y1 = b.y - a.y;
	// Stores coefficient of X
	// direction of vector A[2]A[0]
	// int X2 = (A[2][0] - A[0][0]);
    double x2 = c.x - a.x;
	// Stores coefficient of Y
	// direction of vector A[2]A[0]
	// int Y2 = (A[2][1] - A[0][1]);
    double y2 = c.y - a.y;

    return x1*y2 - y1*x2;
    // Return +ve if notch i.e. reflex angle made by BA to BC in CCW

	// Return cross product
	// return (X1 * Y2 - Y1 * X2);
}

double angle_calcp(pair<double, double> a, pair<double, double> b, pair<double, double> c) {   

    // remove comments
	// Stores coefficient of X
	// direction of vector A[1]A[0]
	// int X1 = (A[1][0] - A[0][0]);
    double x1 = b.first - a.first;
	// Stores coefficient of Y
	// direction of vector A[1]A[0]
	// int Y1 = (A[1][1] - A[0][1]);
    double y1 = b.second - a.second;
	// Stores coefficient of X
	// direction of vector A[2]A[0]
	// int X2 = (A[2][0] - A[0][0]);
    double x2 = c.first - a.first;
	// Stores coefficient of Y
	// direction of vector A[2]A[0]
	// int Y2 = (A[2][1] - A[0][1]);
    double y2 = c.second - a.second;

    return x1*y2 - y1*x2;
    // Return +ve if notch i.e. reflex angle made by BA to BC in CCW

	// Return cross product
	// return (X1 * Y2 - Y1 * X2);
}

DCEL dcelCopier(DCEL d){
    return d;
}

Vertex vCopier(Vertex v){
    return v;
}

DCEL make_DCEL(vector<pair<double, double>> vertex_input){
    DCEL poly;
    int num_v = vertex_input.size();
    
    for(int i = 0; i < num_v; i++){
        // Vertex* temp = (Vertex*) malloc(sizeof(Vertex));
        Vertex *temp = new Vertex;
        temp->key = i;
        temp->x = vertex_input[i].first;
        temp->y = vertex_input[i].second;
        // cout << "hello" << endl;
        poly.vertices.pb(temp);
    }

    for(int i = 0; i < 2 * num_v; i += 2){
        HalfEdge *htemp = new HalfEdge;
        HalfEdge *htwintemp = new HalfEdge;
        htemp->origin = poly.vertices[i/2];
        htemp->origin_vkey = i/2;

        if(i/2 + 1 < num_v){
            htwintemp->origin = poly.vertices[i/2 + 1];
            htwintemp->origin_vkey = i/2 + 1;
        }

        else{
            htwintemp->origin = poly.vertices[0];
            htwintemp->origin_vkey = 0;
        }

        poly.halfEdges.pb(htemp);
        poly.halfEdges.pb(htwintemp);
        // for(auto j : poly.halfEdges){
        //     cout << j.origin->key << " " << j.origin->x << " " << 
        // }
        // cout << "Polyedge Size: " << poly.halfEdges.size() << endl;

        poly.halfEdges.rbegin()[1]->twin = poly.halfEdges.back();
        poly.halfEdges.back()->twin = poly.halfEdges.rbegin()[1];

        poly.vertices[(i/2)]->incidentEdges.pb(poly.halfEdges.rbegin()[1]);
        if((i/2)+1 < num_v) poly.vertices[(i/2) + 1]->incidentEdges.pb(poly.halfEdges.back());
        else poly.vertices[0]->incidentEdges.pb(poly.halfEdges.back());
        
    }

    for(int i = 0; i < 2 * num_v; i += 2) {
        if(i >= 2) {
            poly.halfEdges[i]->prev = poly.halfEdges[i - 2];
            poly.halfEdges[i - 2]->next = poly.halfEdges[i];

            // poly.halfEdges[i + 1]->next = poly.halfEdges[i - 1];
            // poly.halfEdges[i - 1]->prev = poly.halfEdges[i + 1];
        }

        if(i == 0) {
            poly.halfEdges[i]->prev = poly.halfEdges[2 * num_v - 2];
            poly.halfEdges[2 * num_v - 2]->next = poly.halfEdges[i];
        }
    }

    for(auto elt : poly.halfEdges) {
        elt -> twin -> next = elt -> prev -> twin;
        elt -> twin -> prev = elt -> next -> twin;
    }


    swap(poly.vertices[0]->incidentEdges[0], poly.vertices[0]->incidentEdges[1]);
     // swapped so that now 0th element of incident edge is twin


    // cout << poly.vertices[0]->incidentEdges[0]->origin << endl;
    // cout << poly.vertices[0]->incidentEdges[1]->origin << endl;
    // cout << poly.vertices[1]->incidentEdges[0]->origin << endl;
    // cout << poly.vertices[1]->incidentEdges[1]->origin << endl;
    // cout << poly.vertices[2]->incidentEdges[0]->origin << endl;
    // cout << poly.vertices[2]->incidentEdges[1]->origin << endl;
    // cout << poly.vertices[0].incidentEdges.size() << endl;
    // cout << poly.vertices[1].incidentEdges.size() << endl;
    // cout << poly.vertices[2].incidentEdges.size() << endl;
    // cout << poly.halfEdges[4].origin->key;
    // for(auto i : poly.halfEdges){
    //     cout << i.origin_vkey<< endl;
    // }
    // cout << poly.halfEdges.size();
    // for(int i = 0; i < 2 * num_v; i++){
    //     poly.vertices[i/2].incidentEdges.pb()
    // }


    // for(auto &i : poly.halfEdges){
    //     cout << "Address of Node: "<< i << endl;
    //     cout << "Origin_vKey: "<< i->origin_vkey<<endl;
    //     cout << i->origin->x << " " << i->origin->y << endl;
    //     cout << "Twin: " << i->twin << endl;
    //     cout << "Key: " << i->origin->key << endl;
    //     cout << "\n";
    // }

    //  for(auto i : poly.vertices){
    //     cout << i->key << endl;
    //     cout << i->x << " " << i->y << endl;
    //     cout << "IncidentEdge: " << i->incidentEdges[0]->origin->key << " " << i->incidentEdges[0]->origin->x<< " " << i->incidentEdges[0]->origin->y << " "<<i->incidentEdges[0]->origin_vkey << endl;
    //     cout << "IncidentEdge: " << i->incidentEdges[1]->origin->key << " " << i->incidentEdges[1]->origin->x<< " " << i->incidentEdges[1]->origin->y << " "<<i->incidentEdges[1]->origin_vkey << endl;
    //     cout << "IncidentEdges 0  Next: " << i->incidentEdges[0]->next->origin->key << " " << i->incidentEdges[0]->next->origin->x<< " " << i->incidentEdges[0]->next->origin->y << " "<< i->incidentEdges[0]->next->origin_vkey << endl;
    //     cout << "incidentEdges 0 Prev:" << i->incidentEdges[0]->prev->origin->key << " " << i->incidentEdges[0]->prev->origin->x<< " " << i->incidentEdges[0]->prev->origin->y << " "<< i->incidentEdges[0]->prev->origin_vkey << endl;
    //     cout << "incidentEdges 0 Twin: " << i->incidentEdges[0]->twin->origin->key << " " << i->incidentEdges[0]->twin->origin->x<< " " << i->incidentEdges[0]->twin->origin->y << " "<< i->incidentEdges[0]->twin->origin_vkey << endl;

    //      cout << "IncidentEdge Next: " << i->incidentEdges[1]->next->origin->key << " " << i->incidentEdges[1]->next->origin->x<< " " << i->incidentEdges[1]->next->origin->y << " "<< i->incidentEdges[1]->next->origin_vkey << endl;
    //     cout << "incidentEdges Prev:" << i->incidentEdges[1]->prev->origin->key << " " << i->incidentEdges[1]->prev->origin->x<< " " << i->incidentEdges[1]->prev->origin->y << " "<< i->incidentEdges[1]->prev->origin_vkey << endl;
    //     cout << "incidentEdges Twin: " << i->incidentEdges[1]->twin->origin->key << " " << i->incidentEdges[1]->twin->origin->x<< " " << i->incidentEdges[1]->twin->origin->y << " "<< i->incidentEdges[1]->twin->origin_vkey << endl;
    // }
    return poly;
}

vector<Vertex*> findNotches(vector<Vertex*> v) {
    vector<Vertex*> notches;
    for(int i = 0; i < v.size(); i++) {
        if(angle_calc(*v[(i - 1 + v.size()) % v.size()], *v[i], *v[(i+1) % v.size()]) > 0)
            notches.pb(v[i]);
    }
    return notches;
}

int main(){
    freopen("ip1.txt", "r", stdin);

    // Vertices are in CCW order
    // Taking Input
    int num_v;
    int num_e;
    cin >> num_v;
    cin >> num_e;

    vector<pair<double, double>> vertex_input;
    for(int i= 0; i < num_v; i++){
        double x_coor, y_coor;
        cin >> x_coor;
        cin >> y_coor;

        // Vertex* temp = (Vertex*) malloc(sizeof(Vertex));
        // Vertex *temp = new Vertex;
        // temp->key = i;
        // temp->x = x_coor;
        // temp->y = y_coor;
        // cout << "hello" << endl;
        // poly.vertices.pb(temp);
        vertex_input.pb(make_pair(x_coor, y_coor));
    }
    // vertex_input.clear();
    // vertex_input = {{-4, 0}, {0, 4}, {4, 0}, {0, 2}};
    // for(auto elt : notches)
    //     cout << elt.first << " " << elt.second << "\n";


    vector<pair<double, double>> vinp_copy = vertex_input;
    DCEL poly = make_DCEL(vertex_input);
    for(auto i : poly.vertices){
        cout << i->key << endl;
        cout << i->x << " " << i->y << endl;
        cout << "IncidentEdge: " << i->incidentEdges[0]->origin->key << " " << i->incidentEdges[0]->origin->x<< " " << i->incidentEdges[0]->origin->y << " "<<i->incidentEdges[0]->origin_vkey << endl;
        cout << "IncidentEdge: " << i->incidentEdges[1]->origin->key << " " << i->incidentEdges[1]->origin->x<< " " << i->incidentEdges[1]->origin->y << " "<<i->incidentEdges[1]->origin_vkey << endl;
        cout << "IncidentEdges 0  Next: " << i->incidentEdges[0]->next->origin->key << " " << i->incidentEdges[0]->next->origin->x<< " " << i->incidentEdges[0]->next->origin->y << " "<< i->incidentEdges[0]->next->origin_vkey << endl;
        cout << "incidentEdges 0 Prev:" << i->incidentEdges[0]->prev->origin->key << " " << i->incidentEdges[0]->prev->origin->x<< " " << i->incidentEdges[0]->prev->origin->y << " "<< i->incidentEdges[0]->prev->origin_vkey << endl;
        cout << "incidentEdges 0 Twin: " << i->incidentEdges[0]->twin->origin->key << " " << i->incidentEdges[0]->twin->origin->x<< " " << i->incidentEdges[0]->twin->origin->y << " "<< i->incidentEdges[0]->twin->origin_vkey << endl;

         cout << "IncidentEdge Next: " << i->incidentEdges[1]->next->origin->key << " " << i->incidentEdges[1]->next->origin->x<< " " << i->incidentEdges[1]->next->origin->y << " "<< i->incidentEdges[1]->next->origin_vkey << endl;
        cout << "incidentEdges Prev:" << i->incidentEdges[1]->prev->origin->key << " " << i->incidentEdges[1]->prev->origin->x<< " " << i->incidentEdges[1]->prev->origin->y << " "<< i->incidentEdges[1]->prev->origin_vkey << endl;
        cout << "incidentEdges Twin: " << i->incidentEdges[1]->twin->origin->key << " " << i->incidentEdges[1]->twin->origin->x<< " " << i->incidentEdges[1]->twin->origin->y << " "<< i->incidentEdges[1]->twin->origin_vkey << endl;
    }

    DCEL poly_copy;
    // poly_copy = DCEL(poly);
    poly_copy = dcelCopier(poly);

    // for(auto elt : findNotches(poly_copy.vertices))
    //     cout << elt->x << " " << elt->y << "\n";

    // poly2 = poly;

    // cout << poly_copy.vertices[0] << " "<<poly.vertices[0] << endl;
    // poly.vertices[0]->x = 1;
    // poly.vertices[0]->y = 1;
    // cout << poly2.vertices[0]->x << endl;
    // cout << poly2.vertices[0]->y << endl;
    vector<Vertex*> P;

    P = poly_copy.vertices;

    vector<vector<Vertex*>> L;
    vector<Vertex*> vec_temp;
    
    vec_temp.pb(poly_copy.vertices[0]);
    L.pb(vec_temp);
    int m = 1;
    int n = num_v;
    // int i;
    while(n > 3){

        Vertex v0, v1, v_i_plus_one;
        //P = poly_copy.vertices, (P here is supposed to be var you defined in lab)
        for(int i = 0; i < poly_copy.vertices.size(); i++)
            if(P[i] -> x == L[m - 1].back() -> x && P[i] -> y == L[m-1].back() -> y) {
                v0 = *P[i];
                break;
            }

        // v0 = (*L[m-1].back());
        // v2 = (*v1.incidentEdges[1]->next->origin);
        v1 = *poly_copy.vertices[1];
        vector<Vertex*> vv_temp;
        vv_temp.pb(&v0);
        vv_temp.pb(&v1);
        L.pb(vv_temp);
        vv_temp.clear();
        
        // L.pb(vv_temp);
        // vv_temp.clear();
        int i = 1;
        // v_i_plus_one = vCopier(*L[m].back()->incidentEdges[1]->next->origin);
        v_i_plus_one = *poly_copy.vertices[i+1];
        //3.3
        while((angle_calc(*poly_copy.vertices[i-1], *poly_copy.vertices[i], *poly_copy.vertices[i+1]) < 0)&&
              (angle_calc(*poly_copy.vertices[i], *poly_copy.vertices[i+1], v0) < 0) &&
              (angle_calc(*poly_copy.vertices[i+1], v0, v1) < 0) &&
              (L[m].size() < n )
              )
            {       //3.3.1
                    L[m].pb(poly_copy.vertices[i+1]);
                    i = i + 1;
                    v_i_plus_one = *poly_copy.vertices[i+1];
            }

        int min_x = L[m][0]->x, max_x = L[m][0]->x;
        int min_y = L[m][0]->y, max_y = L[m][0]->y;
        for(auto elt : L[m]) {
            if(elt->x < min_x)
                min_x = elt->x;
            if(elt->x > max_x)
                max_x = elt->x;
            if(elt->y < min_y)
                min_y = elt->y;
            if(elt->y > max_y)
                max_y = elt->y;
        }

    // 3.4 
        if(L[m].size() != P.size()) {
            vector<Vertex*> P_minus_Lm;

            vector<Vertex*> notches;
            notches = findNotches(poly_copy.vertices);

            // vinp_copy defined above, will be changed as iters progress
            //3.4.1
            vector<Vertex*> LPVS;
            for(int i = 1; i < L[m].size() - 1; i++) {
                for(auto elt : poly_copy.vertices) {
                    if(elt->x == L[m][i]->x && elt->y == L[m][i]->y)
                        continue;
                    else{
                        for(auto notch : notches) {
                            if(elt->x == notch->x && elt->y == notch->y){
                                LPVS.pb(elt);
                                break;
                            }
                        }
                    }
                }
            }

            //3.4.2
            while(LPVS.size() > 0) {
                
                bool backward = false;
                while(!backward && LPVS.size()){
                    // Resume from here
                    Vertex* v = new Vertex;
                    v = LPVS[0];

                }
            }

        }
    
        if(vertexEqual(*L[m].back(), v1)){
            

        }





    }

   

    return 0;
}



// checking if a point is inside a polygon or not  https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
// sort poly vertices everytime on basis of key
// superscript indexing on vertices indexes inside L[m]
// L[m] should have pointer from the immediate previous copy of poly
// dcel next and prev how, when not simple polygon