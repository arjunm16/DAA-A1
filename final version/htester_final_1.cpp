#include <bits/stdc++.h>
#include "dcel_hvg.hpp"
#define pb push_back

using namespace std;

bool vertexEqual(Vertex v1, Vertex v2){
    if(v1.x == v2.x && v1.y == v2.y)
        return true;

    return false;
}

struct Point{
    double x, y;
};

struct line{
    Point p1, p2;
};

void printLm(vector<vector<Vertex*>> L){
    int count = 0;
    for(auto i : L){
        cout << "count : " << count << endl;
        for(auto j : i){
            cout << j->x << " " << j->y << endl;
        }
        count++;
    }
}

tuple<double, double, double> givecoef(double x1, double y1, double x2, double y2){
    
	double A, B, C;
	if(x2 == x1){
		A = 1;
		B = 0;
		C = -x1;
		tuple<double, double, double> tup(A, B, C);
		return tup;
	}
    A = (y1 - y2) / (x2 - x1);
	B = 1;
    C = (x1 * (y2 - y1)/(x2 - x1)) - y1;

  	tuple<double, double, double> tup(A, B, C);
	return tup;
}


bool cmp(pair<double, double> &a, pair<double, double> &b)
{

    if (a.first == b.first){
        int ans = b.second > a.second;
        return ans;
    }
    int ans = b.first > a.first;
    return ans;
}

int cw(pair<double, double> &a, pair<double, double> &b, pair<double, double> &c) {

    double p = a.first * (b.second - c.second)+ b.first * (c.second - a.second)+ c.first * (a.second - b.second);
    int ans = p < 0;
    return ans;
}

int ccw(pair<double, double>&a, pair<double, double>&b, pair<double, double>&c){
    double p = a.first * (b.second - c.second)+ b.first * (c.second - a.second)+ c.first * (a.second - b.second);
    int ans = p > 0;
    return ans;
}

void iterover(vector<pair<double, double>> &v, int &n, pair<double, double> &p1, pair<double, double> &p2, vector<pair<double, double>> &t1, vector<pair<double, double>> &t2){
    for (int i = 1; i < n; i++){
         if (i == n - 1 || !ccw(p1, v[i], p2)) {
            while (t1.size() > 1 && ccw(t1[t1.size() - 2],t1[t1.size() - 1],v[i])) {
                t1.pop_back();
            }
            t1.push_back(v[i]);
        }

        if (i == n - 1 || !cw(p1, v[i], p2)) {
            while (t2.size() > 1 && cw(t2[t2.size() - 2],t2[t2.size() - 1],v[i])) {
                t2.pop_back();
            }
            t2.push_back(v[i]);
        }
    }
}

vector<pair<double, double>> convexHull(vector<pair<double, double>> &v)
{
    sort(v.begin(), v.end(), cmp);
    int n = v.size();
    if (n < 4) return v;
    //Taking first point
    pair<double, double> p1 = v[0];
    // Taking last point
    pair<double, double> p2 = v[n - 1];

    vector<pair<double, double>> u, d;

    // Insert StartingEnding Points
    u.push_back(p1);
    d.push_back(p1);

    // Iterate over points

    iterover(v, n, p1, p2, u, d);

    // Combine upper and lower half
    int i = d.size() - 2;
    while(i){
        u.push_back(d[i]);
        i--;
    }
    // Remove duplicate points
    u.resize(unique(u.begin(),u.end())- u.begin());

    // Return the points on Convex Hull
    return u;
}


bool isInside(vector<pair<double, double> > points,
            pair<double, double> query)
{
    points.push_back(query);
    points = convexHull(points);
    for (auto x : points) {
        if (x == query)
            return 0;
    }
    return 1;
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

DCEL make_DCEL(vector<pair<double, double>> vertex_input){
    DCEL poly;
    int num_v = vertex_input.size();
    
    for(int i = 0; i < num_v; i++){
        Vertex *temp = new Vertex;
        temp->key = i;
        temp->x = vertex_input[i].first;
        temp->y = vertex_input[i].second;
        poly.vertices.pb(temp);
    }

    for(int i = 0; i<num_v; i++)
        poly.vertices[i]->next = poly.vertices[(i + 1) % num_v]; 
    
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
        poly.halfEdges.rbegin()[1]->twin = poly.halfEdges.back();
        poly.halfEdges.back()->twin = poly.halfEdges.rbegin()[1];

        poly.vertices[(i/2)]->incidentEdges.pb(poly.halfEdges.rbegin()[1]);
        if((i/2) + 1 < num_v) poly.vertices[(i/2) + 1]->incidentEdges.pb(poly.halfEdges.back());
        else poly.vertices[0]->incidentEdges.pb(poly.halfEdges.back());
        
    }

    for(int i = 0; i < 2 * num_v; i += 2) {
        if(i >= 2) {
            poly.halfEdges[i]->prev = poly.halfEdges[i - 2];
            poly.halfEdges[i - 2]->next = poly.halfEdges[i];
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

    return poly;
}

vector<Vertex*> findNotches(vector<Vertex*> v) {
    vector<Vertex*> notches;
    for(int i = 0; i < v.size(); i++) {
        if(angle_calc(*v[(i - 1 + v.size()) % v.size()], *v[i], *v[(i + 1) % v.size()]) > 0)
            notches.pb(v[i]);
    }
    return notches;
}

void addDiag(DCEL &d, Vertex* v1, Vertex* v2){
    // cout << "\nmaking diagonal from (" << v1->x << ", " << v1->y << ") to (" << v2->x << ", " << v2->y << ")\n";
    // cout << "i.e. " << v1->key << " to " << v2->key << endl;
    for(auto i : d.vertices){
        if(i->x == v1->x && i->y == v1->y){
            HalfEdge* h1 = new HalfEdge;
            h1->origin = v1;       
            h1->origin_vkey = v1->key;
            h1->prev = v1->incidentEdges[1]->twin;
            v1->incidentEdges[1]->twin->next = h1;
            h1->next = v2->incidentEdges[0];
            v2->incidentEdges[0]->prev = h1;

            HalfEdge* h2 = new HalfEdge; // twin of h1
            h2->origin = v2;
            h2->origin_vkey = v2->key;
            h2->prev = v2->incidentEdges[1]->twin;
            v2->incidentEdges[1]->twin->next = h2;
            h2->next = v1->incidentEdges[0];
            v1->incidentEdges[0]->prev = h2;

            h1->twin = h2;
            h2->twin = h1;

            
            d.diags.pb(h2);
            d.diags.pb(h1);
            break;
        }
    }
}

int main(int argc, char *argv[]) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ifstream infile;
    infile.open(argv[1]);
    
    // freopen(argv[1], "r", stdin);

    // Vertices are in CCW order
    // Taking Input
    int num_v;
    infile >> num_v;

    vector<pair<double, double>> vertex_input;
    for(int i = 0; i < num_v; i++){
        double x_coor, y_coor;
        infile >> x_coor;
        infile >> y_coor;
        vertex_input.pb(make_pair(x_coor, y_coor));
    }

    infile.close();

    vector<pair<double, double>> vinp_copy = vertex_input;

    DCEL poly = make_DCEL(vertex_input);

    cout << findNotches(poly.vertices).size() << endl;

    DCEL poly_copy;
    poly_copy = dcelCopier(poly);

    vector<vector<Vertex*>> L;
    vector<Vertex*> vec_temp;
    vector<int> convexIndex;

    vec_temp.pb(poly_copy.vertices[0]); 
    L.pb(vec_temp);
    
    int m = 1;
    int n = num_v;
    
    while(n > 3){
        // cout << "Iteration : " << m << endl;
        vector<Vertex*> P;
        P = poly_copy.vertices;
        Vertex *v0, *v1, *v_i_plus_one;
        //search in polygon for v0 and then put those in v0 and v1
        for(auto i : P) {
            if((i->x == L[m-1].back()->x) && (i->y == L[m-1].back()->y)) {
                v0 = i;
                break;
            }
        }
        // cout << v0->x << " " << v0->y << endl;

        v1 = v0->next;
        // cout << v1->x << " " << v1->y << endl;

        vector<Vertex*> vv_temp;
        vv_temp.pb(v0);
        vv_temp.pb(v1);
        L.pb(vv_temp);
        vv_temp.clear();

        int iter = 1;
        v_i_plus_one = v1->next;
        while((angle_calc(*L[m][iter-1], *L[m][iter], *v_i_plus_one) <= 0)&&
              (angle_calc(*L[m][iter], *v_i_plus_one, *v0) <= 0) &&
              (angle_calc(*v_i_plus_one, *v0, *v1) <= 0) &&
              (L[m].size() < n )
              )
            {       //3.3.1
                    L[m].pb(v_i_plus_one);
                    iter = iter + 1;
                    v_i_plus_one = v_i_plus_one->next;
                    // cout << "v_i_plus_one "<<v_i_plus_one->x << " " << v_i_plus_one->y << endl;

            }
        // cout << "L[m] size after angle while: " << L[m].size() << endl;
        // cout << "L[m] : \n";
        // for(auto i : L[m]){
        //     cout << i->x << " " << i->y << endl;
        // }
        // cout << "L[m] finished\n";
        
    // 3.4 
        if(L[m].size() != P.size()) {
            vector<Vertex*> notches;
            notches = findNotches(P);
            // vinp_copy defined above, will be changed as iters progress
            //3.4.1
            vector<Vertex*> LPVS;
            // making LPVS vector
            for(auto elt : P) {
                // cout << elt->x << " " << elt->y << endl;
                bool flag = false;
                for(int k = 0; k < L[m].size(); k++){
                    if(elt->x == L[m][k]->x && elt->y == L[m][k]->y) {
                        flag = true;
                        break;
                    }
                } 
                // cout << "flag : " << flag << endl;
                if(!flag) {
                    // if notch add else do nothing
                    for(auto notch : notches)
                        if(elt->x == notch->x && elt->y == notch->y) {
                            // cout << "pushed\n";
                            LPVS.pb(elt);
                            break;
                        }
                }
            }

            //obtain smallest rect R with sides || to axes and containing all vertices of Lm
            //3.4.2
            // cout << "LPVS size: " << LPVS.size() << endl;
            // for(auto i : LPVS){
            //     cout << i->x << " " << i->y << endl;
            // }
            // int LPVS_counter = 1;
            while(LPVS.size() > 0) {
                // cout << "LPVS_counter : " << LPVS_counter << endl;
                bool backward = false;
                while(!backward && (LPVS.size() > 0)){
                    while(LPVS.size() != 0){
                        Vertex* v = new Vertex;
                        v = LPVS[0];
                        double min_x = 100000.0, max_x = -100000.0;
                        double min_y = 100000.0, max_y = -100000.0;
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
                        // cout << "min_x max_x min_y max_y : \n";
                        // cout << min_x << " " <<  max_x << " " <<  min_y << " " <<  max_y << " " << endl;
                        // cout << "v : " << endl;
                        // cout << v->x << " "<<v->y << endl;
                        if(v->x >= min_x && v->x <= max_x && v->y >= min_y && v->y <= max_y){
                            // cout << "breaked in rectangle: \n";
                            break;
                        }
                        else if(LPVS.size() == 0){
                            break;
                        }
                        else{
                            LPVS.erase(LPVS.begin());
                        }
                    }
                    if(LPVS.size() > 0){
                        vector<pair<double, double>> pts;
                        for(auto i : L[m]){
                            pts.pb(make_pair(i->x, i->y));
                        }
                        bool inside = isInside(pts, make_pair(LPVS[0]->x, LPVS[0]->y));
                        // cout << "inside: " << inside << endl;
                        if(inside){
                            vector<Vertex*> VTR;
                            tuple<double, double, double> coef = givecoef(L[m][0]->x, L[m][0]->y, LPVS[0]->x, LPVS[0]->y);
                            double sign = (get<0>(coef) * L[m].back()->x) + (get<1>(coef) * L[m].back()->y) + (get<2>(coef));
                            for(auto elt: L[m]){
                                if(sign * ((get<0>(coef) * elt->x) + (get<1>(coef) * elt->y) + (get<2>(coef))) > 0){
                                    VTR.pb(elt);
                                }
                            }
                            // cout << "VTR: \n";
                            // for(auto i: VTR){
                            //     cout << i->x << " " << i->y << endl;
                            // }
                            vector<Vertex*> Lm_VTR;
                            // cout << "size Lm_VTR: " << Lm_VTR.size() << endl; 
                            // printLm(L);
                            for(auto elt : L[m]){
                                bool flag1 = false;
                                for(auto j : VTR){
                                    if(elt->x == j->x && elt->y == j->y){
                                        flag1 = true;
                                        break;
                                    }    
                                }
                                if(!flag1){
                                    // cout << "pushed\n";
                                    // cout << elt->x << " " << elt->y << endl;
                                    Lm_VTR.pb(elt);
                                }
                            }
                            L.pop_back();
                            L.push_back(Lm_VTR);
                            // cout << "Lm_VTR: \n";
                            // for(auto i : Lm_VTR){
                            //     cout << i->x << " " << i->y << endl;
                            // }
                            backward = true;
                            LPVS.erase(LPVS.begin()); 
                        }
                        else{
                            LPVS.erase(LPVS.begin()); 
                        }
                    }
                }
                // LPVS_counter++;
            }
        }

        // cout << "vertexEqual: "<<vertexEqual(*L[m].back(), *v1) << endl;
        if(!vertexEqual(*L[m].back(), *v1)){
            convexIndex.pb(m);
            
            Vertex *src = NULL, *dest = NULL;
            for(auto vert : poly.vertices)
                if(vert->x == L[m][0]->x && vert->y == L[m][0]->y) {
                    src = vert;
                    break;
                }
            
            for(auto vert : poly.vertices)
                if(vert->x == L[m].back()->x && vert->y == L[m].back()->y) {
                    dest = vert;
                    break;
                }
            // cout << "source " << src->x << " " << src->y << endl;
            // cout << "dest " << dest->x << " " << dest->y << endl;
            // addDiag(poly, *L[m][0], *(L[m].back()));
            bool dExists = false;

            for(auto elt : poly.diags) {
                if(src->x == elt->origin->x && src->y == elt->origin->y && dest->x == elt->twin->origin->x && dest->y == elt->twin->origin->y) {
                    dExists = true;
                    break;
                }
            }

            if(!dExists)
                addDiag(poly, src, dest);

            vector<Vertex*> updatedPolyVerts;
            for(auto elt : P){
                bool flag2 = false;
                for(int j = 1; j < L[m].size() - 1; j++){
                    if(elt->x == L[m][j]->x && elt->y == L[m][j]->y){
                        flag2 = true;
                        break;
                    }    
                }
                if(!flag2)
                    updatedPolyVerts.pb(elt);
            }

            // cout << "updatedPolyVerts : " << endl;
            // for(auto i : updatedPolyVerts){
            //     cout << i->x << " " << i->y << endl;
            // }
            vector<pair<double, double>> vinpnew;
            for(auto elt : updatedPolyVerts)
                vinpnew.pb(make_pair(elt->x, elt->y));
            
            poly_copy = make_DCEL(vinpnew);
            n = n - L[m].size() + 2;
        }
        m = m + 1;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << endl; //prints us

    // string fname = "op";

    ofstream outfile;
    outfile.open(("op" + string(argv[1])).c_str());
    // freopen((fname + string(argv[1])).c_str(), "w", stdout);

    for(int i = 0; i < convexIndex.size(); i++){
        for(auto j : L[convexIndex[i]]){
            outfile << j->x << " " << j->y << " ";
        }
        outfile << endl;
    }

    for(auto l: poly.vertices)
        outfile << l->x << " " << l->y << " ";
    outfile << endl;

    outfile.close();

    // fname = "diag";
    outfile.open(("diag" + string(argv[1])).c_str());

    int c = 0;
    for(auto i : poly.diags) {
        outfile << ++c << endl;
        outfile << "origin: key: " << i->origin->key << " -> (" << i->origin->x << ", " << i->origin->y << ")" << endl;
        outfile << "next origin: key: " << i->next->origin->key << " -> (" << i->next->origin->x << ", " << i->next->origin->y << ")" << endl;
        outfile << "prev origin: key: " << i->prev->origin->key << " -> (" << i->prev->origin->x << ", " << i->prev->origin->y << ")" << endl;
        outfile << "twin origin: key: " << i->twin->origin->key << " -> (" << i->twin->origin->x << ", " << i->twin->origin->y << ")" << endl;
        outfile << "twin next origin: key: " << i->twin->next->origin->key << " -> (" << i->twin->next->origin->x << ", " << i->twin->next->origin->y << ")" << endl;
        outfile << "twin prev origin: key: " << i->twin->prev->origin->key << " -> (" << i->twin->prev->origin->x << ", " << i->twin->prev->origin->y << ")" << endl;
        outfile << endl;
    }

    outfile << "Number of elements in poly.diags: " << poly.diags.size() << endl;
    outfile << "Number of diagonals: " << poly.diags.size() / 2;
    return 0;
}