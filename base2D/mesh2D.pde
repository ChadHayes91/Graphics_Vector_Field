// PLANAR TRIANGLE MESH
// Jarek Rossignac, Nov 6, 2019

class MESH {
  
    // VERTICES
    int nv=0, maxnv = 1000, numNewVerts=0;  
    pt[] G = new pt [maxnv];               // location of vertex
    vec[] F = new vec [maxnv];             // vector at vertex
    pt[] newVerts = new pt [maxnv];        // new vertices to add for retriangulation
    vec[] newVecs = new vec [maxnv];       // the vectors at each of the new vertices
    pt[] backupG = new pt [maxnv];         // stores the vertices before retriangulation
    vec[] backupF = new vec [maxnv];       // stores the vectors before retriangulation
    int backupNV = 0;                      // the number of vertices before retriangulation
    
    // TRIANGLES 
    int nt = 0, maxnt = maxnv*2;                           
    boolean[] isInterior = new boolean[maxnv];   // true if the triangle at the specified index is not a border triangle
    boolean[] touched_tris = new boolean[maxnv]; // true if the trangle at the specified index has been touched by the current trace
    int[] isFollowedBy = new int [maxnv];        // the index of triangle that follows the triangle at the specified index on the computed path
    int mouseT = -1;                             // the triangle that the mouse was last in when a mouse trace was requested
    
    // CORNERS 
    int c=0;    // current corner                                                              
    int nc = 0; // corner count
    int[] V = new int [3*maxnt];   // Corner table c.v
    int[] O = new int [3*maxnt];   // Corner table o.v
    
    // TRACES
    pts[] traces = new pts[maxnv];             // Point cload objects for each trace
    int[] traceStartCorners = new int [maxnv]; // Start locations for each trace
    int numTraces = 0;
    pts[] bgTraces = new pts[maxnv]; // All backgroundTraces
    int numBgTraces = 0;
    
    // DRAWING
    float traceWeight = 10;
    float stepSize = 3;
     
  MESH() {for (int i=0; i<maxnv; i++) {G[i]=P(); F[i]=V(); newVerts[i] = P();}}; // declare all points and vectors
  void reset() {
    nv=0; 
    nt=0; 
    nc=0; 
    for (int i=0; i<maxnv; i++) {
      touched_tris[i] = false;
    }
  }     // removes all vertices and triangles
  void resetToBackups() { G = backupG; F = backupF; nv = backupNV; numNewVerts=0; }
  void loadVertices(pt[] P, int n) {nv=0; for (int i=0; i<n; i++) addVertex(P[i]);}
  void writeVerticesTo(pts P) {for (int i=0; i<nv; i++) P.G[i].setTo(G[i]);}
  void addVertex(pt P) { G[nv++].setTo(P); }                                             // adds a vertex to vertex table G
  void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt=nc/3; }     // adds triangle (i,j,k) to V table
  void addVertexPVfromPP(pt A, pt B) {G[nv].setTo(A); F[nv++].setTo(V(A,B));}                                             // adds a vertex to vertex table G
  void loadFromPTS(pts P) {int n=P.nv; nv=0; for (int i=0; i<n; i+=2) addVertexPVfromPP(P.G[i],P.G[i+1]);}

  // CORNER OPERATORS
  int t (int c) {int r=int(c/3); return(r);}                   // triangle of corner c
  int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);}         // next corner
  int p (int c) {int r=3*int(c/3)+(c+2)%3; return(r);}         // previous corner
  int v (int c) {return V[c];}                                // vertex of c
  int o (int c) {return O[c];}                                // opposite corner
  int l (int c) {return o(n(c));}                             // left
  int s (int c) {return n(o(n(c)));}                             // left
  int u (int c) {return p(o(p(c)));}                             // left
  int r (int c) {return o(p(c));}                             // right
  pt g (int c) {return G[V[c]];}                             // shortcut to get the point where the vertex v(c) of corner c is located
  vec f (int c) {return F[V[c]];}                             // shortcut to get the vector of the vertex v(c) of corner c 
  pt cg(int c) {return P(0.8,g(c),0.1,g(p(c)),0.1,g(n(c)));}   // computes offset location of point at corner c

  boolean nb(int c) {return(O[c]!=c);};  // not a border corner
  boolean bord(int c) {return(O[c]==c);};  // a border corner
  int firstBorderCorner() {int i=0; while(nb(i) && i<nc) i++; return i;}
  pt firstBorderEdgeMidPoint() {int fbc = M.firstBorderCorner(); return P(g(p(fbc)),g(n(fbc)));}
  
  /* Returns the following integer array:
   *   [0] = number of triangles touched
   *   [1] = length of the trace in number of steps
   */
  int[] tracePathFromMidEdgeFacingCorner(int sc, boolean draw, boolean save_touched, boolean save_new_verts) // sc = start corner
  {
    pt P = P(g(p(sc)),g(n(sc))); // start at midpoint of edge facing sc
    int c = sc;
    pt Q = P();
    pen(brown,traceWeight); noFill(); 
    int[] traceInfo = traceFromPoint(P, t(c), draw, save_touched, save_new_verts, false, 0); // trace forwards
    int[] traceInfo2 = traceFromPoint(P, t(c), draw, save_touched, save_new_verts, false, 0, true); // trace backwards
    if (traceInfo[2] != traceInfo2[2] && save_new_verts && traceInfo2[0] > 0) isFollowedBy[traceInfo2[2]] = traceInfo[2];
    P.setTo(Q);
    return traceInfo;
  }
  
  /* Call this version for saving new traces
   *   Returns the following integer array:
   *   [0] = number of triangles touched
   *   [1] = length of the trace in number of steps
   */
  int[] tracePathFromMidEdgeFacingCorner(int sc, boolean draw, boolean save_touched, boolean save_new_verts, boolean save_trace, int index) // sc = start corner
  {
    pt P = P(g(p(sc)),g(n(sc))); // start at midpoint of edge facing sc
    int c = sc;
    pt Q = P();
    int[] traceInfo = traceFromPoint(P, t(c), draw, save_touched, save_new_verts, save_trace, index);
    int[] traceInfo2 = traceFromPoint(P, t(c), draw, save_touched, save_new_verts, save_trace, index, true);
    if (traceInfo[2] != traceInfo2[2] && save_new_verts && traceInfo2[0] > 0) isFollowedBy[traceInfo2[2]] = traceInfo[2];
    P.setTo(Q);
    return traceInfo;
  }
  
  /* Returns the following integer array:
   *   [0] = number of triangles touched
   *   [1] = length of the trace in number of steps
   *   [2] = the first triangle touched
   */
  int[] traceFromPoint(pt P, int t, boolean draw, boolean save_touched, boolean save_new_verts, boolean save_trace, int index)
  {
    int iters_in_tri = 0;
    int max_iters = 2000;
    int hard_count = 0;
    int hard_cap = 5000;
    int first_t = t;
    int last_t = t;
    int numTouched = 0;
    pt entryPt = P;
    boolean[] backup = touched_tris.clone();
    while (getRawContainingTriangle(P) != -1 && iters_in_tri < max_iters && (hard_count < hard_cap)) {
      if (hard_count > 0) touched_tris[t] = true;
      vec P_ = getVectorAt(P, t);
      pt q = P(P, P_);
      vec q_ = getVectorAt(q, getContainingTriangle(q, t));
      pt r = P(q, M(q_));
      vec S = V(P, r);
      pt p_next = P(q, 0.5, S);
      vec dir = V(P, p_next).normalize();
      p_next = P(P, stepSize, dir);
      if (draw || tracingMouse) {
        line(P.x, P.y, p_next.x, p_next.y);
      }
      if (save_trace) traces[index].addPt(P);
      P = p_next;
      t = getContainingTriangle(P, t);      
      if (hard_count == 0) first_t = t;
      hard_count++;
      if (((last_t != t && hard_count == 1) || last_t == t) && getRawContainingTriangle(P) != -1) {
        if (last_t != t && hard_count == 1 && touched_tris[t]) break;
        iters_in_tri++;
        last_t = t;
      }
      else {
        if (save_new_verts) {
          pt mid = getTraceMidpoint(entryPt, iters_in_tri / 2, last_t, false);
          newVerts[last_t] = mid;
          newVecs[last_t] = getVectorAt(P, last_t);
          if (t != -1 && !touched_tris[t]) isFollowedBy[last_t] = t;
          else isFollowedBy[last_t] = -1;
          numNewVerts++;
        }
        if (save_trace) traces[index].addPt(P);
        iters_in_tri = 0;
        if (touched_tris[t]) break;
        last_t = t;
        entryPt = P;
        numTouched++;
      }
    }
    if (iters_in_tri == max_iters) {
      if (save_new_verts) {
        newVerts[t] = P;
        newVecs[t] = getVectorAt(P, t);
        isFollowedBy[t] = -1;
        numNewVerts++;
      }
    }
    if (!save_touched) touched_tris = backup.clone();
    int[] traceInfo = new int[3];
    traceInfo[0] = numTouched;
    traceInfo[1] = hard_count;
    traceInfo[2] = first_t;
    return traceInfo;
  }
  
  /* Backwards traces
   * Returns the following integer array:
   *   [0] = number of triangles touched
   *   [1] = length of the trace in number of steps
   *   [2] = the first triangle touched
   */
  int[] traceFromPoint(pt P, int t, boolean draw, boolean save_touched, boolean save_new_verts, boolean save_trace, int index, boolean backwards)
  {
    int original_t = t;
    int iters_in_tri = 0;
    int max_iters = 2000;
    int hard_count = 0;
    int hard_cap = 5000;
    int last_t = t;
    int first_t = t;
    int numTouched = 0;
    pt entryPt = P;
    boolean[] backup = touched_tris.clone();
    while (getRawContainingTriangle(P) != -1 && iters_in_tri < max_iters && (hard_count < hard_cap)) {
      if (hard_count > 0) touched_tris[t] = true;
      vec P_ = W(-1, getVectorAt(P, t));
      pt q = P(P, P_);
      vec q_ = getVectorAt(q, getContainingTriangle(q, t));
      pt r = P(q, M(q_));
      vec S = V(P, r);
      pt p_next = P(q, 0.5, S);
      vec dir = V(P, p_next).normalize();
      p_next = P(P, stepSize, dir);
      if (draw || tracingMouse) {
        line(P.x, P.y, p_next.x, p_next.y);
      }
      if (save_trace) traces[index].prependPt(P);
      P = p_next;
      t = getContainingTriangle(P, t);
      if (hard_count == 0) first_t = t;
      hard_count++;
      if (((last_t != t && hard_count == 1) || last_t == t) && getRawContainingTriangle(P) != -1) {
        if (last_t != t && hard_count == 1 && touched_tris[t]) break;
        iters_in_tri++;
        last_t = t;
      }
      else {
        if (save_new_verts && t != original_t) {
          pt mid = getTraceMidpoint(entryPt, iters_in_tri / 2, last_t, backwards);
          newVerts[last_t] = mid;
          newVecs[last_t] = getVectorAt(P, last_t);
          if (t != -1 && !touched_tris[t]) isFollowedBy[t] = last_t;
          numNewVerts++;
        }
        if (save_trace) traces[index].prependPt(P);
        iters_in_tri = 0;
        if (touched_tris[t]) break;
        last_t = t;
        entryPt = P;
        numTouched++;
      }
    }
    if (iters_in_tri == max_iters) {
      if (save_new_verts) {
        newVerts[t] = P;
        newVecs[t] = getVectorAt(P, t);
        isFollowedBy[t] = -1;
        numNewVerts++;
      }
    }
    if (!save_touched) touched_tris = backup.clone();
    if (save_trace) traces[index].removeFirstPt();
    int[] traceInfo = new int[3];
    traceInfo[0] = numTouched;
    traceInfo[1] = hard_count;
    traceInfo[2] = first_t;
    return traceInfo;
  }
  
  // Traces a given number of iterations along a trace, returns the end point.  Used to get the midpoint of a trace in a triangle.
  pt getTraceMidpoint(pt start, int iters, int t, boolean backwards) {
    int i = 0;
    pt P = P(start);
    while (i < iters) {
      vec P_ = getVectorAt(P, t);
      if (backwards) P_ = W(-1, P_);
      pt q = P(P, P_);
      vec q_ = getVectorAt(q, t);
      pt r = P(q, M(q_));
      vec S = V(P, r);
      pt p_next = P(q, 0.5, S);
      vec dir = V(P, p_next).normalize();
      p_next = P(P, stepSize, dir);
      P = p_next;
      i++;
    }
    return P;
  }
  
  // Finds the best locations for traces to start.  Saves them to traceStartCorners.
  void calculateTraceStarts() {
    // Reset trace data structures
    numTraces = 0;
    
    // Loop until all triangles have been reached
    while (!allTrisTouched()) {
      int max_touched = -1; // max number of triangles touched by any one trace
      int longestTrace = -1; // length of the longest trace 
      int bestCorner = -1; // Start corner for longest trace
      
      // loop over all triangles
      for (int t = 0; t < nt; t++) {
        if (touched_tris[t] == false) {
          for (int i = 0; i < 3; i++) {
            int[] traceInfo = tracePathFromMidEdgeFacingCorner((t * 3) + i, false, false, false);
            if (traceInfo[0] > max_touched) {
              max_touched = traceInfo[0];
              bestCorner = t * 3 + i;
              longestTrace = traceInfo[1];
            } else if (traceInfo[0] == max_touched) {
              if (traceInfo[1] > longestTrace) {
                max_touched = traceInfo[0];
                bestCorner = t * 3 + i;
                longestTrace = traceInfo[1];
              }
            }
          }
        }
      }
      tracePathFromMidEdgeFacingCorner(bestCorner, false, true, false); // Updates touched_tris
      traceStartCorners[numTraces] = bestCorner; // Saves the start corner for later drawing
      numTraces++;
    }
    
    // Reset touched_tris to false
    for (int i = 0; i < nt; i++) {
      touched_tris[i] = false;
    }
  }
  
  // Traces from the given point a set number of iterations. //<>//
  void traceIterations(pt P, int t, int iters, boolean save, boolean backwards) {
    if (save) {
      bgTraces[numBgTraces] = new pts();
      bgTraces[numBgTraces].declare();
      bgTraces[numBgTraces].empty();
    }
    for (int i = 0; i < iters; i++) {
      if (save) bgTraces[numBgTraces].addPt(P);
      vec P_ = getVectorAt(P, t);
      if (backwards) P_ = W(-1, P_);
      pt q = P(P, P_);
      vec q_ = getVectorAt(q, getContainingTriangle(q, t));
      pt r = P(q, M(q_));
      vec S = V(P, r);
      pt p_next = P(q, 0.5, S);
      vec dir = V(P, p_next).normalize();
      p_next = P(P, stepSize, dir);
      P = p_next;
      t = getContainingTriangle(P, t);
      if (getRawContainingTriangle(P) == -1) break;
    }
    if (save) numBgTraces++;
  }
  
  // Uses a naive approach to set trace start locations.
  void calculateNaiveTraceStarts() {    
    // loop over all triangles
    for (int t = 0; t < nt; t++) {
      if (touched_tris[t] == false) {
        int longest = -1;
        int bestCorner = 0;
        for (int i = 0; i < 3; i++) {
          int[] traceInfo = tracePathFromMidEdgeFacingCorner((t * 3) + i, false, false, false);
          if (traceInfo[1] > longest) {
            longest = traceInfo[1];
            bestCorner = t * 3 + i;
          }
        }
        tracePathFromMidEdgeFacingCorner(bestCorner, false, true, false); // Updates touched_tris
        traceStartCorners[numTraces] = bestCorner; // Saves the start corner for later drawing
        numTraces++;
      }
    }
    
    // Reset touched_tris to false
    for (int i = 0; i < nt; i++) {
      touched_tris[i] = false;
    }
  }
  
  // Saves the traces from traceStartCorners to the traces array as point clouds to be displayed later.
  void saveTraces() {
    for (int i = 0; i < numTraces; i++) {
      traces[i] = new pts();
      traces[i].declare();
      tracePathFromMidEdgeFacingCorner(traceStartCorners[i], false, true, false, true, i);
    }
    
    // Reset touched_tris to false
    for (int i = 0; i < nt; i++) {
      touched_tris[i] = false;
    }
  }
  
  // Calculates new vertices to add to the triangle mesh.  These are at the midpoint of each trace in a triangle.
  void calculateNewVerts() {
    // Reset newVerts data structure
    numNewVerts = 0;
    
    for (int i = 0; i < numTraces; i++) {
      tracePathFromMidEdgeFacingCorner(traceStartCorners[i], false, true, true);
    }
    
    // Reset touched_tris to false
    for (int i = 0; i < nt; i++) {
      touched_tris[i] = false;
    }
  }
  
  void drawTraceStarts() {
    for (int i = 0; i < numTraces; i++) {
      tracePathFromMidEdgeFacingCorner(traceStartCorners[i], true, false, false);
    }
  }
  
  // Checks if all triangles in the mesh have been touched by the current set of traces.
  boolean allTrisTouched()
  {
    for (int i = 0; i < nt; i++) {
      if (!touched_tris[i]) return false;
    }
    return true;
  }
  
  // Returns the vector at any given location in the mesh.  Needs a triangle as input to know what parent vectors to use.
  vec getVectorAt(pt P, int t)
  {
    int[] verts = getTrianglePoints(t);
    int A = verts[0];
    int B = verts[1];
    int C = verts[2];
    vec AP = V(G[A], P);
    vec AC = V(G[A], G[C]);
    vec AB = V(G[A], G[B]);
    float b = det(AP, AC) / det(AB, AC);
    float c = det(AP, AB) / det(AC, AB);
    float a = 1 - b - c;
    vec temp = W(a, F[A], b, F[B]);
    return W(temp, c, F[C]);
  }
  
  // Return the index of the triangle that contains the given point.  If no triangle contains the point, it returns the value
  // passed to it in the t parameter.
  int getContainingTriangle(pt P, int t) {
    int tri = getRawContainingTriangle(P);
    if (tri == -1) return t;
    else return tri;
  }
  
  // Return the index of the triangle that contains the given point.  If no triangle contains the point, returns -1.
  int getRawContainingTriangle(pt P) {
    for (int i = 0; i < nt; i++) {
      int[] verts = getTrianglePoints(i);
      if (ptInTriangle(P, G[verts[0]], G[verts[1]], G[verts[2]])) return i;
    }
    return -1;
  }
  
  // Returns the ids of the three vertices of triangle t.
  int[] getTrianglePoints(int t) {
    int[] pts = new int[3];
    int c = t*3;
    pts[0] = v(c);
    pts[1] = v(c+1);
    pts[2] = v(c+2);
    return pts;
  }
  
  // RETRIANGULATION
  void retriangulate() {
    // Add vertices
    int old_nt = nt;
    int old_nv = nv;
    reset();
    for (int i = 0; i < old_nv; i++) {
      addVertex(backupG[i]);
      F[nv-1] = backupF[i];
      
    }
    for (int i = 0; i < old_nt; i++) {
      addVertex(newVerts[i]);
      F[nv-1] = newVecs[i];
    }
    
    // Triangulate
    triangulate();
    isRetriangulated = true;
  }
  
  // Sets up the backup data for reinstantiating the old triangluation if needed.
  void saveData() {
    // Copy G and F into backups
    for (int i = 0; i < maxnv; i++) {
      backupG[i] = P(G[i]);
      backupF[i] = V(F[i]);
    }
    backupNV = nv;
  }
  
  // Checks if the given point is in the newVerts array.
  boolean ptInNewVerts(pt P) {
    for (int i = 0; i < numNewVerts; i++) {
      if (P.x == newVerts[i].x && P.y == newVerts[i].y) {
        return true;
      }
    }
    return false;
  }
  
  // Checks if the two points given are consecutive newVerts on the same trace.
  boolean ptsFollow(pt A, pt B) {
    if (ptInNewVerts(A) && ptInNewVerts(B)) {
      // Find which triangle A and B were in
      int a = -1;
      int b = -1;
      for (int i = 0; i < numNewVerts; i++) {
        if (A.x == newVerts[i].x && A.y == newVerts[i].y) {
          a = i;
        }
        if (B.x == newVerts[i].x && B.y == newVerts[i].y) {
          b = i;
        }
      }
      return isFollowedBy[a] == b || isFollowedBy[b] == a;
    }
    else return false;
  }
  
  // CURRENT CORNER OPERATORS
  void next() {c=n(c);}
  void previous() {c=p(c);}
  void opposite() {c=o(c);}
  void left() {c=l(c);}
  void right() {c=r(c);}
  void swing() {c=s(c);} 
  void unswing() {c=u(c);} 
  void printCorner() {println("c = "+c);}
  
  // DISPLAY
  void showCurrentCorner(float r) { show(cg(c),r); };   // renders corner c 
  void showEdge(int c) {edge( g(p(c)),g(n(c))); };  // draws edge of t(c) opposite to corner c
  void showVertices(float r) {for (int v=0; v<nv; v++) show(G[v],r); }                          // shows all vertices 
  void showBorderVertices(float r) {for (int v=0; v<nv; v++) if(!isInterior[v]) show(G[v],r);} // shows only border vertices              
  void showInteriorVertices(float r) {for (int v=0; v<nv; v++) if(isInterior[v]) show(G[v],r); }   // shows interior vertices 
  // draws all triangles (edges, or filled)
  void showTriangles() { 
    for (int c=0; c<nc; c+=3) {
      noFill();
      pen(black,2);
      if (ptsFollow(g(c), g(c+1))) pen(black, 7);
      line(g(c).x, g(c).y, g(c+1).x, g(c+1).y);
      pen(black, 2);
      if (ptsFollow(g(c+1), g(c+2))) pen(black, 7);
      line(g(c+1).x, g(c+1).y, g(c+2).x, g(c+2).y);
      pen(black, 2);
      if (ptsFollow(g(c), g(c+2))) pen(black, 7);
      line(g(c).x, g(c).y, g(c+2).x, g(c+2).y);
    }
  } 
  void showEdges() {for (int i=0; i<nc; i++) showEdge(i); };         // draws all edges of mesh twice
  void showBorderEdges() {for (int i=0; i<nc; i++) {if (bord(i)) {showEdge(i);}; }; };         // draws all border edges of mesh
  void showNonBorderEdges() {for (int i=0; i<nc; i++) {if (!bord(i)) {showEdge(i);}; }; };         // draws all border edges of mesh
  void showVerticesAndVectors(float r) {for (int v=0; v<nv; v++) {show(G[v],r); arrow(G[v],F[v]);}}   // shows all vertices 
  void drawArrows() 
    {
    stroke(blue); 
    for (int v=0; v<nv; v++) 
      { 
      fill(blue); arrow(G[v],F[v]);
      fill(white); 
      show(G[v],13); 
      fill(black); 
      if(v<10) label(G[v],str(v));  
      else label(G[v],V(-1,0),str(v)); 
      }
    noFill();
    }
  void showCorner(int c, float r) { if(bord(c)) show(cg(c),1.5*r); else show(cg(c),r); };   // renders corner c 
  void showCorners(float r) 
    {
    noStroke(); 
    for (int c=0; c<nc; c+=3) 
      {
      fill(red); showCorner(c,r); 
      fill(dgreen); showCorner(c+1,r); 
      fill(blue); showCorner(c+2,r);
      } 
    }
  void fillTriangle(int t) {
    int[] verts = getTrianglePoints(t);
    fill(yellow);
    beginShape();
    vertex(G[verts[0]].x, G[verts[0]].y);
    vertex(G[verts[1]].x, G[verts[1]].y);
    vertex(G[verts[2]].x, G[verts[2]].y);
    endShape();
  }
  void showCoPVectors() {
    for (int c = 0; c < nc; c+=3) {
      pt center = P(g(c), g(c+1), g(c+2));
      vec v = getVectorAt(center, int(c/3));
      fill(black);
      stroke(black);
      show(center, 3);
      show(Arrow(center, v));
    }
  }
  void setMouseTriangle() {
    pt P = P(mouseX, mouseY);
    int t = getRawContainingTriangle(P);
    mouseT = t;
  }
  void traceMouseTriangle() {
    if (mouseT != -1) {
      int c = mouseT*3;
      tracingMouse = true;
      traceFromPoint(P(g(c), g(c+1), g(c+2)), mouseT, false, false, false, false, 0);
      tracingMouse = false;
    }
  }
  void showTraces() {
    pen(brown, traceWeight);
    for (int i = 0; i < numTraces; i++) {
      traces[i].drawOpenCurve();
    }
    
    // Reset touched_tris
    for (int i = 0; i < nt; i++) {
      touched_tris[i] = false;
    }
  }
  void showNewVerts() {
    for (int i = 0; i < numNewVerts; i++) {
      fill(red);
      show(newVerts[i], 5);
    }
  }
  
  // Shows a bunch of faint traces throughout the field
  void showBackgroundTraces() {
    for (int i = 0; i < numBgTraces; i++) {
      bgTraces[i].drawOpenCurve();
    }
  }
  
  void saveBgTraces() {
    float inc = 20; // Number of pixels between each trace //<>//
    int len = 500; // The max length of each trace
    numBgTraces = 0;
    for (int c = 0; c < nc; c++) {
      if (bord(c)) {
        pt p1 = g(n(c));
        pt p2 = g(p(c));
        vec dir = V(p1, p2).normalize();
        float d = 0;
        while (d < d(p1, p2)) {
          pt origin = P(p1, d, dir);
          pen(grey, 1);
          traceIterations(origin, t(c), len, true, false);
          d += inc;
        }
      }
    }
  }

  // DISPLAY
  void classifyVertices() 
    { 
    for (int v=0; v<nv; v++) isInterior[v]=true;
    for (int c=0; c<nc; c++) if(bord(c)) isInterior[v(n(c))]=false;
    }               

 void triangulate() {     // performs Delaunay triangulation using a quartic algorithm
   c=0;                   // to reset current corner
   pt X = new pt(0,0);
   float r=1;
   for (int i=0; i<nv-2; i++) for (int j=i+1; j<nv-1; j++) for (int k=j+1; k<nv; k++) {
     X=CircumCenter(G[i],G[j],G[k]);  r = d(X,G[i]);
     boolean found=false; 
     for (int m=0; m<nv; m++) if ((m!=i)&&(m!=j)&&(m!=k)&&(d(X,G[m])<=r)) found=true;  
     if (!found) {
       if (cw(G[i],G[j],G[k])) addTriangle(i,j,k); 
       else addTriangle(i,k,j);
     };
   }; 
 }  

  void computeO() {   // slow method to set the O table from the V table, assumes consistent orientation of tirangles
    for (int i=0; i<3*nt; i++) {O[i]=i;};  // init O table to -1: has no opposite (i.e. is a border corner)
    for (int i=0; i<3*nt; i++) {  for (int j=i+1; j<3*nt; j++) {       // for each corner i, for each other corner j
      if( (v(n(i))==v(p(j))) && (v(p(i))==v(n(j))) ) {O[i]=j; O[j]=i;};};}; // make i and j opposite if they match 
   }
  
  void computeOfast() // faster method for computing O
    {                                          
    int nIC [] = new int [maxnv];                            // number of incident corners on each vertex
    println("COMPUTING O: nv="+nv +", nt="+nt +", nc="+nc );
    int maxValence=0;
    for (int c=0; c<nc; c++) {O[c]=c;};                      // init O table to -1: has no opposite (i.e. is a border corner)
    for (int v=0; v<nv; v++) {nIC[v]=0; };                    // init the valence value for each vertex to 0
    for (int c=0; c<nc; c++) {nIC[v(c)]++;}                   // computes vertex valences
    for (int v=0; v<nv; v++) {if(nIC[v]>maxValence) {maxValence=nIC[v]; };};  println(" Max valence = "+maxValence+". "); // computes and prints maximum valence 
    int IC [][] = new int [maxnv][maxValence];                 // declares 2D table to hold incident corners (htis can be folded into a 1D table !!!!!)
    for (int v=0; v<nv; v++) {nIC[v]=0; };                     // resets the valence of each vertex to 0 . It will be sued as a counter of incident corners.
    for (int c=0; c<nc; c++) {IC[v(c)][nIC[v(c)]++]=c;}        // appends incident corners to corresponding vertices     
    for (int c=0; c<nc; c++) {                                 // for each corner c
      for (int i=0; i<nIC[v(p(c))]; i++) {                     // for each incident corner a of the vertex of the previous corner of c
        int a = IC[v(p(c))][i];      
        for (int j=0; j<nIC[v(n(c))]; j++) {                   // for each other corner b in the list of incident corners to the previous corner of c
           int b = IC[v(n(c))][j];
           if ((b==n(a))&&(c!=n(b))) {O[c]=n(b); O[n(b)]=c; };  // if a and b have matching opposite edges, make them opposite
           };
        };
      };
    } // end computeO
  
  pt triCenter(int c) {return P(g(c),g(n(c)),g(p(c))); }  // returns center of mass of triangle of corner c
  pt triCircumcenter(int c) {return CircumCenter(g(c),g(n(c)),g(p(c))); }  // returns circumcenter of triangle of corner c

  void smoothenInterior() { // even interior vertiex locations
    pt[] Gn = new pt[nv];
    int[] sum = new int[nv];
    for (int v=0; v<nv; v++) sum[v]=0;
    for (int v=0; v<nv; v++) Gn[v]=P(0,0);
    for (int c=0; c<3*nt; c++) 
      {
      float d=d(g(n(c)),g(p(c))); 
      Gn[v(c)].add(d,P(g(n(c)),g(p(c)))); 
      sum[v(c)]+=d;
      }
    for (int v=0; v<nv; v++) Gn[v].scale(1./sum[v]);
    for (int v=0; v<nv; v++) if(isInterior[v]) G[v].translateTowards(.1,Gn[v]);
    }

  } // end of MESH
