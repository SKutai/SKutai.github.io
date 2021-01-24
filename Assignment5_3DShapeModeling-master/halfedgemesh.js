/**
 * Skeleton code implementation for a half-edge mesh
 */

let vec3 = glMatrix.vec3;

class HEdge {
    /**
     * Constructor for a half-edge
     */
    constructor() {
        this.head = null; // Head vertex (Type HVertex)
        this.face = null; // Left face (Type HFace)
        this.pair = null; // Half edge on opposite face (Type HEdge)
        this.prev = null; // Previous half edge in CCW order around left face (Type HEdge)
        this.next = null; // Next half edge in CCW order around left face (Type HEdge)

        this.checked = false;
    }

    /**
     * Return a list of the two vertices attached to this edge,
     * or an empty list if one of them has not yet been initialized
     * 
     * @returns {list} A 2-element list of HVertex objects corresponding
     *                  to vertices attached to this edge
     */
    getVertices() {
        let ret = [];
        if (!(this.head === null) && !(this.prev === null)) {
            if (!(this.prev.head === null)) {
                ret = [this.head, this.prev.head];
            }
        }
        return ret;
    }
}

class HFace {
    /**
     * Constructor for a half-edge face
     */
    constructor() {
        this.h = null; // Any HEdge on this face (Type HEdge)

        this.checked = false;
    }
    
    /**
     * Get a list of half-edges involved with this face
     * 
     * @returns {list} A list of HEdge objects corresponding
     *                 to edges at the boundary of this face
     */
    getEdges() {
        if (this.h === null) {
            return [];
        }
        // TODO: Fill this in
        let edge0 = this.h;
        let nextEdge = this.h;
        let edges = [];
        
        do{
            edges.push(nextEdge);
            nextEdge = nextEdge.next;
        }while(nextEdge != edge0);

        return edges;
    }

    /**
     * Get a list of vertices attached to this face
     * 
     * @returns {list} A list of HVertex objects corresponding
     *                 to vertices on this face
     */
    getVertices() {
        if (this.h === null) {
            return [];
        }
        let h = this.h.next;
        let vertices = [this.h.head];
        while (h != this.h) {
            vertices.push(h.head);
            h = h.next;
        }
        return vertices;
    }

    /**
     * Compute the area of this face
     * 
     * @returns {float} The area of this face
     */
    getArea() {
        let area = 0.0;
        // TODO: Fill this in (you can use mini assignment 1 to help)
        // Remember, there are n-2 triangles in an n-gon
        let edge0 = this.h;
        let vecU = vec3.create();
        let vecV = vec3.create();
        let crossProd = vec3.create();

        let a = this.h.head.pos;
        this.h = this.h.next;

        let b = this.h.head.pos;
        this.h = this.h.next;

        let c = this.h.head.pos;

        vec3.sub(vecU, b, a);
        vec3.sub(vecV, c, a);

        do{
            vec3.cross(crossProd, vecU, vecV);
            area += (vec3.len(crossProd) / 2);
            b = c;
            this.h = this.h.next;
            c = this.h.head.pos;
            vec3.sub(vecU, b, a);
            vec3.sub(vecV, c, a);

        }while(this.h != edge0);
        return area;
    }

    /**
     * Get the normal of this face, assuming it is flat
     * 
     * @returns {vec3} The normal of this face
     */
    getNormal() {
        let normal = vec3.create();
        // TODO: Fill this in
        let edge0 = this.h;
        let vecU = vec3.create();
        let vecV = vec3.create();
        let crossProd = vec3.create();

        let a = this.h.head.pos;
        this.h = this.h.next;

        let b = this.h.head.pos;
        this.h = this.h.next;

        let c = this.h.head.pos;

        vec3.sub(vecU, b, a);
        vec3.sub(vecV, c, a);

        vec3.cross(crossProd, vecU, vecV);
        vec3.normalize(normal, crossProd);

        this.h = edge0;

        return normal;
    }
}

class HVertex {
    /**
     * Constructor for a half-edge vertex
     * @param {glMatrix.vec3} pos Position of this vertex
     * @param {glMatrix.vec3} color Color of this vertex
     */
    constructor(pos, color) {
        this.pos = pos; // Position of this vertex (Type vec3)
        this.color = color; // Color of this vertex (Type vec3)
        this.h = null; // Any hedge on this vertex (Type Hedge)

        this.checked = false;
    }

    /**
     * Compute the vertices that are attached to this
     * vertex by an edge
     * 
     * @returns {list} List of HVertex objects corresponding
     *                 to the attached vertices
     */
    getVertexNeighbors() {
        if (this.h === null) {
            return [];
        }
        // TODO: Fill this in (a do while loop works well)
        let neighbor = this.h.head;
        let neighbors = [];
        let edge0 = this.h;

        do{
            neighbors.push(neighbor);
            this.h = this.h.pair.next;
            neighbor = this.h.head;
        }while(this.h != edge0);

        return neighbors;
    }

    /**
     * Compute the faces of which this vertex is a member
     * 
     * @returns {list} A list of HFace objects corresponding
     *                  to the incident faces
     */
    getAttachedFaces() {
        if (this.h === null) {
            return [];
        }
        // TODO: Fill this in (this is very similar to getVertexNeighbors)
        let faces = [];
        let edge0 = this.h;

        do{
            faces.push(this.h.face);
            this.h = this.h.pair.next;
        }while(this.h != edge0);

        return faces;
    }

    /**
     * Compute the normal of this vertex as an area-weighted
     * average of the normals of the faces attached to this vertex
     * 
     * @returns {vec3} The estimated normal
     */
    getNormal() {
        // TODO: This is a dummy value that makes all normals (1, 0, 0)
        // just so something shows up.  You should actually initialize it to (0, 0, 0)
        // and accumulate neighboring face normals in a loop for your average
        let normal = vec3.fromValues(0, 0, 0); 
        // TODO: Fill this in 
        // Hint: use this.getAttachedFaces(), face.getArea(), and face.getNormal() to help
        let out = vec3.create();
        let faces = this.getAttachedFaces();

        for(let i = 0; i < faces.length; i++){
            if(faces[i] != null){
                normal = vec3.scaleAndAdd(out, normal, faces[i].getNormal(), faces[i].getArea());
            }
        }

        normal = vec3.normalize(out, normal);

        return normal;
    }
}


class HedgeMesh extends PolyMesh {
    /**
     * @returns {I} A NumTrisx3 Uint16Array of indices into the vertex array
     */
    getTriangleIndices() {
        let NumTris = 0;
        let allvs = [];
        for (let i = 0; i < this.faces.length; i++) {
            let vsi = this.faces[i].getVertices();
            allvs.push(vsi.map(function(v){
                return v.ID;
            }));
            NumTris += vsi.length - 2;
        }
        let I = new Uint16Array(NumTris*3);
        let i = 0;
        let faceIdx = 0;
        //Now copy over the triangle indices
        while (i < NumTris) {
            let verts = allvs[faceIdx]
            for (let t = 0; t < verts.length - 2; t++) {
                I[i*3] = verts[0];
                I[i*3+1] = verts[t+1];
                I[i*3+2] = verts[t+2];
                i++;
            }
            faceIdx++;
        }
        return I;
    }

    /**
     * @returns {I} A NEdgesx2 Uint16Array of indices into the vertex array
     */
    getEdgeIndices() {
        let I = [];
        for (let i = 0; i < this.edges.length; i++) {
            let vs = this.edges[i].getVertices();
            for (let k = 0; k < vs.length; k++) {
                I.push(vs[k].ID);
            }
        }
        return new Uint16Array(I);
    }

    /**
     * Given two vertex objects representing an edge,
     * and a face to the left of that edge, initialize
     * a half edge object and add it to the list of edges
     * 
     * @param {HVertex} v1 First vertex on edge
     * @param {HVertex} v2 Second vertex on edge
     * @param {HFace} face Face to the left of edge
     * 
     * @returns {HEdge} The constructed half edge
     */
    addHalfEdge(v1, v2, face) {
        const hedge = new HEdge();
        hedge.head = v2; // Points to head vertex of edge
        hedge.face = face;
        v1.h = hedge; // Let tail vertex point to this edge
        this.edges.push(hedge);
        return hedge;
    }

    /////////////////////////////////////////////////////////////
    ////                INPUT/OUTPUT METHODS                /////
    /////////////////////////////////////////////////////////////

    /**
     * Load in an OFF file from lines and convert into
     * half edge mesh format. Crucially, this function assumes
     * a consistently oriented mesh with vertices specified 
     * in CCW order
     */
    loadFileFromLines(lines) {
        // Step 1: Consistently orient faces using
        // the basic mesh structure and copy over the result
        const origMesh = new BasicMesh();
        origMesh.loadFileFromLines(lines);
        origMesh.consistentlyOrientFaces();
        origMesh.subtractCentroid();
        const res = {'vertices':[], 'colors':[], 'faces':[]};
        for (let i = 0; i < origMesh.vertices.length; i++) {
            res['vertices'].push(origMesh.vertices[i].pos);
            res['colors'].push(origMesh.vertices[i].color);
        }
        for (let i = 0; i < origMesh.faces.length; i++) {
            // These faces should now be consistently oriented
            const vs = origMesh.faces[i].getVertices();
            res['faces'].push(vs.map(
                function(v) {
                    return v.ID;
                }
            ));
        }

        // Step 1.5: Clear previous mesh
        this.vertices.length = 0;
        this.edges.length = 0;
        this.faces.length = 0;

        // Step 2: Add vertices
        for (let i = 0; i < res['vertices'].length; i++) {
            let V = new HVertex(res['vertices'][i], res['colors'][i]);
            V.ID = this.vertices.length;
            this.vertices.push(V);
        }

        let str2Hedge = {};
        // Step 3: Add faces and halfedges
        for (let i = 0; i < res['faces'].length; i++) {
            const face = new HFace();
            this.faces.push(face);
            let vertsi = [];
            for (let k = 0; k < res['faces'][i].length; k++) {
                vertsi.push(this.vertices[res['faces'][i][k]]);
            }

            // Add halfedges
            for (let k = 0; k < vertsi.length; k++) {
                const v1 = vertsi[k];
                const v2 = vertsi[(k+1)%vertsi.length];
                // Add each half edge
                const hedge = this.addHalfEdge(v1, v2, face);
                // Store half edge in hash table
                let key = v1.ID+"_"+v2.ID;
                str2Hedge[key] = hedge;
                face.h = hedge;
            }

            // Link edges together around face in CCW order
            // assuming each vertex points to the half edge
            // starting at that vertex
            // (which addHalfEdge has just done)
            for (let k = 0; k < vertsi.length; k++) {
                vertsi[k].h.next = vertsi[(k+1)%vertsi.length].h;
                vertsi[(k+1)%vertsi.length].h.prev = vertsi[k].h;
            }
        }

        // Step 4: Add links between opposite half edges if 
        // they exist.  Otherwise, it is a boundary edge, so
        // add a half edge with a null face on the other side
        let boundaryEdges = {}; // Index boundary edges by their tail
        for (const key in str2Hedge) {
            const v1v2 = key.split("_");
            let h1 = str2Hedge[key];
            const other = v1v2[1]+"_"+v1v2[0];
            if (other in str2Hedge) {
                h1.pair = str2Hedge[other];
            }
            else {
                let h2 = new HEdge();
                h1.pair = h2;
                h2.pair = h1;
                h2.head = this.vertices[v1v2[0]];
                boundaryEdges[v1v2[1]] = h2;
                this.edges.push(h2);
            }
        }

        // Step 5: Link boundary edges
        for (let key in boundaryEdges) {
            let e = boundaryEdges[key];
            if (e.next === null) {
                let e2 = boundaryEdges[e.head.ID];
                e.next = e2;
                e2.prev = e;
            }
        }

        console.log("Initialized half edge mesh with " + 
                    this.vertices.length + " vertices, " + 
                    this.edges.length + " half edges, " + 
                    this.faces.length + " faces");

        this.needsDisplayUpdate = true;
    }


    /////////////////////////////////////////////////////////////
    ////                  GEOMETRIC TASKS                   /////
    /////////////////////////////////////////////////////////////
    
    /**
     * Move each vertex along its normal by a factor
     * 
     * @param {float} fac Move each vertex position by this
     *                    factor of its normal.
     *                    If positive, the mesh will inflate.
     *                    If negative, the mesh will deflate.
     */
    inflateDeflate(fac) {
        // Loop through all the vertices in the mesh
        for (let vertex of this.vertices) {
            // TODO: Fill this in
            vec3.scaleAndAdd(vertex.pos, vertex.pos, vertex.getNormal(), fac);
        }
        this.needsDisplayUpdate = true;
    }

    /**
     * Compute the mean vector from all of this vertex's neighbors
     * to the vertex.  If smoothing, subtract this vector off.
     * If sharpening, add this vector on
     * 
     * @param {boolean} smooth If true, smooth.  If false, sharpen
     */
    laplacianSmoothSharpen(smooth) {
        // TODO: Fill this in
        let laplacianVectors = [];
        let neighbors = [];
        let neighborCount = 0;
        let mainVertPos = vec3.create();
        let laplacianVector = vec3.create();
        let edge = vec3.create();

        for (let vertex of this.vertices) {
            mainVertPos = vertex.pos;
            neighbors = vertex.getVertexNeighbors();
            neighborCount = neighbors.length;
            laplacianVector = vec3.fromValues(0,0,0);
            
            for (let nearVertex of neighbors){
                vec3.sub(edge, nearVertex.pos, mainVertPos);
                vec3.add(laplacianVector, laplacianVector, edge);
            }

            vec3.scale(laplacianVector, laplacianVector, (1/neighborCount));
            laplacianVectors.push(laplacianVector);
        }

        for(let i = 0; i < this.vertices.length; i++){
            if(smooth){
                vec3.add(this.vertices[i].pos, this.vertices[i].pos, laplacianVectors[i]);
            }
            else{
                vec3.scale(laplacianVectors[i], laplacianVectors[i], -1);
                vec3.add(this.vertices[i].pos, this.vertices[i].pos, laplacianVectors[i]);
            }
        }

        this.needsDisplayUpdate = true;
    }

    /** Apply some creative per-vertex warp */
    warp() {
        // TODO: Fill this in

        this.needsDisplayUpdate = true;
    }


    /////////////////////////////////////////////////////////////
    ////                  TOPOLOGICAL TASKS                 /////
    /////////////////////////////////////////////////////////////
    /**
     * Return a list of boundary cycles
     * 
     * @returns {list} A list of cycles, each of which is
     *                 its own list of HEdge objects corresponding
     *                 to a unique cycle
     */
    getBoundaryCycles() {
        let cycles = [];
        // TODO: Fill this in (hint: Add a variable to an edge which
        // stores whether this edge has been checked yet)
        let faceless = [];
        let edge0 = this.edges[0];

        for (let halfEdge of this.edges) {
            if(halfEdge.checked == false){
                if(halfEdge.face === null){
                    edge0 = halfEdge;
                    do{
                        halfEdge.checked = true;
                        faceless.push(halfEdge);
                        halfEdge = halfEdge.next;
                    }while(halfEdge != edge0);
                    cycles.push(faceless);
                }
                halfEdge.checked = true;
            }
        }

        for (let halfEdge of this.edges) {
            halfEdge.checked = false;
        }

        return cycles;
    }

    /**
     * Compute the genus of this mesh if it is watertight.
     * If it is not watertight, return -1
     * 
     * @returns {int} genus if watertight, or -1 if not
     */
    getGenus() {
        let genus = -1;
        // TODO: Fill this in (hint: there are two half edges for every edge!)

        if(this.getBoundaryCycles().length == 0){
            let Nv = this.vertices.length;
            let Ne = this.edges.length / 2;
            let Nf = this.faces.length;
            
            let X = Nv - Ne + Nf;

            genus = Math.floor((-X + 2) / 2);
            
        }

        return genus;

    }

    /**
     * Fill in the boundary cycles with triangles.  The mesh
     * should be watertight at the end
     */
    fillHoles() {
        // TODO: Fill this in

        this.needsDisplayUpdate = true;
    }



    /////////////////////////////////////////////////////////////
    ////                MESH CREATION TASKS                 /////
    /////////////////////////////////////////////////////////////
    
    /**
     * Truncate the mesh by slicing off the tips of each vertex
     * @param {float} fac The amount to go down each edge from the vertex
     *                    (should be between 0 and 1)
     */
    truncate(fac) {
        // TODO: Fill this in

        this.needsDisplayUpdate = true;
    }

    /**
     * Perform a linear subdivision of the mesh
     */
    subdivideLinear() {
        // TODO: Fill this in

        this.needsDisplayUpdate = true;
    }

    /** 
     * Perform Loop subdivision on the mesh
     */
    subdivideLoop() {
        // TODO: Fill this in

        this.needsDisplayUpdate = true;
    }
}