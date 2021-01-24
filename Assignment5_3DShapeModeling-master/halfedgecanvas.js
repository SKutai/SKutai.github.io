DRAW_OFFSET = 1.001; // Draw points at a slight multiplicative offset for debugging
LINE_WIDTH = 4;

class HalfEdgeCanvas extends BaseCanvas {
    /**
     * @param {DOM Element} glcanvas Handle to HTML where the glcanvas resides
     * @param {string} shadersrelpath Path to the folder that contains the shaders,
     *                                relative to where the constructor is being called
     * @param {antialias} boolean Whether antialiasing is enabled (true by default)
     */
    constructor(glcanvas, shadersrelpath, antialias) {
        super(glcanvas, shadersrelpath, antialias);
        this.mesh = new HedgeMesh();
        this.camera = new MousePolarCamera(glcanvas.width, glcanvas.height);
        // Setup drawer object for debugging.  It is undefined until
        // the pointColorShader is ready
        let canvas = this;
        if (!('shaderReady' in this.shaders.pointColorShader)) {
            this.shaders.pointColorShader.then(function() {
                canvas.drawer = new SimpleDrawer(canvas.gl, canvas.shaders.pointColorShader);
            })
        }
        else {
            this.drawer = new SimpleDrawer(this.gl, this.shaders.pointColorShader);
        }
        this.faceDrawer = new BasicMesh();
        this.setupMenus();
    }

    centerCamera() {
        this.camera.centerOnMesh(this.mesh);
    }

    /**
     * Setup GUI menus for students to test their functions
     */
    setupMenus() {
        this.gui = new dat.GUI();
        const gui = this.gui;
        let canvas = this;
        let simpleRepaint = function() {
            requestAnimFrame(canvas.repaint.bind(canvas));
        }
        // Mesh display options menu
        this.drawEdges = false;
        this.drawNormals = false;
        this.drawVertices = false;
        let meshOpts = gui.addFolder('Mesh Display Options');
        ['drawEdges', 'drawNormals', 'drawPoints'].forEach(
            function(s) {
                let evt = meshOpts.add(canvas, s);
                evt.onChange(simpleRepaint);
            }
        );
    
        // Add tests for mesh traversal
        let travTests = gui.addFolder('Traversal Tests');
        let vertTests = travTests.addFolder('Vertex Tests');
        this.vertTest = {'vertexIndex':0, 'showResult':false, 'type':'getVertexNeighbors'};
        this.vertTest.gui = vertTests;
        vertTests.add(this.vertTest, 'vertexIndex').onChange(function(i) {
            let maxi = canvas.mesh.vertices.length-1;
            if (i >= maxi) {
                canvas.vertTest.vertexIndex = maxi;
            }
            if (i < 0) {
                canvas.vertTest.vertexIndex = 0;
            }
            simpleRepaint();
        });
        vertTests.add(this.vertTest, 'type', ['getVertexNeighbors', 'getAttachedFaces']).onChange(simpleRepaint);
        vertTests.add(this.vertTest, 'showResult').onChange(simpleRepaint);
    
        let faceTests = travTests.addFolder('Face Tests');
        this.faceTest = {'faceIndex':0, 'showResult':false, 'type':'getEdges'};
        this.faceTest.gui = faceTests;
        faceTests.add(this.faceTest, 'faceIndex').onChange(function(i) {
            let maxi = canvas.mesh.faces.length-1;
            if (i >= maxi) {
                canvas.faceTest.faceIndex = maxi;
            }
            if (i < 0) {
                canvas.faceTest.faceIndex = 0;
            }
            simpleRepaint();
        });
        faceTests.add(this.faceTest, 'type', ['getEdges']).onChange(simpleRepaint);
        faceTests.add(this.faceTest, 'showResult').onChange(simpleRepaint);
    
        let geomMenu = gui.addFolder("Geometric Tasks");
        this.inflationFac = 0.1;
        this.inflateDeflate = function() {
            canvas.mesh.inflateDeflate(canvas.inflationFac);
            simpleRepaint();
        }
        this.laplacianSmooth = function() {
            canvas.mesh.laplacianSmoothSharpen(true);
            simpleRepaint();
        }
        this.laplacianSharpen = function() {
            canvas.mesh.laplacianSmoothSharpen(false);
            simpleRepaint();
        }
        this.warp = function() {
            canvas.mesh.warp();
            simpleRepaint();
        }
        geomMenu.add(this, 'inflationFac', -1, 1);
        geomMenu.add(this, 'inflateDeflate');
        geomMenu.add(this, 'laplacianSmooth');
        geomMenu.add(this, 'laplacianSharpen');
        geomMenu.add(this, 'warp');
    
        let topoMenu = gui.addFolder("Topological Tasks");
        this.showBoundaries = false;
        this.genus = -1;
        topoMenu.add(this, 'showBoundaries').onChange(simpleRepaint);
        topoMenu.add(this, 'genus').listen();
        this.fillHoles = function() {
            canvas.mesh.fillHoles();
            simpleRepaint();
        }
        topoMenu.add(this, 'fillHoles');
    
        let creationMenu = gui.addFolder("Mesh Creation");
        this.truncationFac = 0.5;
        creationMenu.add(this, 'truncationFac', 0.01, 0.99);
        this.truncate = function() {
            canvas.mesh.truncate(this.truncationFac);
            simpleRepaint();
        }
        creationMenu.add(this, 'truncate');
    
        this.subdivideLinear = function() {
            canvas.mesh.subdivideLinear();
            simpleRepaint();
        }
        creationMenu.add(this, 'subdivideLinear');
        this.subdivideLoop = function() {
            canvas.mesh.subdivideLoop();
            simpleRepaint();
        }
        creationMenu.add(this, 'subdivideLoop');
    
        gui.add(this.mesh, 'saveOffFile').onChange(simpleRepaint);
        simpleRepaint();
    }

    /**
     * Code for showing result of face-related traversals
     * @param {HEdge} edge Half-edge
     * @param {glMatrix.vec3} color color
     */
    drawEdge(edge, color) {
        let drawer = this.drawer;
        if (color === undefined) {
            color = [1, 0, 1];
        }
        let vs = edge.getVertices();
        let offset = DRAW_OFFSET;
        for (let k = 0; k < 2; k++) {
            if (k == 0) {
                offset = DRAW_OFFSET;
            }
            else {
                offset = 1.0/DRAW_OFFSET;
            }
            let p1 = glMatrix.vec3.create();
            let p2 = glMatrix.vec3.create();
            glMatrix.vec3.scale(p1, vs[0].pos, DRAW_OFFSET);
            glMatrix.vec3.scale(p2, vs[1].pos, DRAW_OFFSET);
            drawer.drawLine(p1, p2, color);
        }
    }

    /**
     * Code for showing result of vertex-related traversals
     */
    drawVertexTraversals() {
        let drawer = this.drawer;
        if (this.vertTest.showResult) {
            let v = this.mesh.vertices[this.vertTest.vertexIndex];
            let p1 = glMatrix.vec3.create();
            let p1_2 = glMatrix.vec3.create();
            glMatrix.vec3.scale(p1, v.pos, DRAW_OFFSET);
            glMatrix.vec3.scale(p1_2, v.pos, 1.0/DRAW_OFFSET);
            drawer.drawPoint(p1, [1, 0, 0]);
            if (this.vertTest.type == 'getVertexNeighbors' || this.vertTest.type == 'getAttachedFaces') {
                let vs = v.getVertexNeighbors();
                for (let i = 0; i < vs.length; i++) {
                    let p2 = glMatrix.vec3.create();
                    let p2_2 = glMatrix.vec3.create();
                    glMatrix.vec3.scale(p2, vs[i].pos, DRAW_OFFSET);
                    glMatrix.vec3.scale(p2_2, vs[i].pos, 1.0/DRAW_OFFSET);
                    drawer.drawLine(p1, p2, [0, 1, 1]);
                    drawer.drawPoint(p2, [0, 1, 1]);
                    drawer.drawLine(p1_2, p2_2, [0, 1, 1]);
                    drawer.drawPoint(p2_2, [0, 1, 1]);
                }
            }
            if (this.vertTest.type == 'getAttachedFaces') {
                // Draw all attached faces in red
                let faces = v.getAttachedFaces();
                this.faceDrawer.clear();
                let idx2V = {}; // Convert from mesh index to new vertex in face drawer
                for (let i = 0; i < faces.length; i++) {
                    const f = faces[i];
                    const vs = f.getVertices();
                    let face = [];
                    let face_2 = [];
                    for (let j = 0; j < vs.length; j++) {
                        const v = vs[j];
                        if (!(v.ID in idx2V)) {
                            let p = glMatrix.vec3.create();
                            let p_2 = glMatrix.vec3.create();
                            glMatrix.vec3.scale(p, v.pos, DRAW_OFFSET);
                            glMatrix.vec3.scale(p_2, v.pos, 1.0/DRAW_OFFSET);
                            const vnew = this.faceDrawer.addVertex(p);
                            const vnew_2 = this.faceDrawer.addVertex(p_2);
                            idx2V[v.ID] = [vnew, vnew_2];
                        }
                        face.push(idx2V[v.ID][0]);
                        face_2.push(idx2V[v.ID][1])
                    }
                    this.faceDrawer.addFace(face);
                    this.faceDrawer.addFace(face_2);
                }
                let materialBefore = this.material;
                this.material = {'ka':[1, 1, 0], 'kd':[1, 1, 0]};
                this.faceDrawer.render(this);
                this.material = materialBefore;
            }
        }
    }

    /**
     * Code for drawing the mesh and debugging information
     */
    repaint() {
        let canvas = this;
        let drawer = this.drawer;
        if (!('shaderReady' in this.shaders.blinnPhong)) {
            // Wait until the promise has resolved, then draw again
            this.shaders.blinnPhong.then(canvas.repaint.bind(canvas));
            return;
        }
        if (!('shaderReady' in this.shaders.pointColorShader)) {
            // Wait until the promise has resolved, then draw again
            this.shaders.pointColorShader.then(canvas.repaint.bind(canvas));
            return;
        }
        let gl = this.gl;
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.lineWidth(LINE_WIDTH);

        this.lights = [{pos:this.camera.pos, color:[1, 1, 1], atten:[1, 0, 0]}];
        this.shaderToUse = this.shaders.blinnPhong;
        this.mesh.render(this);

        
        drawer.reset();
        this.drawVertexTraversals();

        if (this.faceTest.showResult) {
            let face = this.mesh.faces[this.faceTest.faceIndex];
            let edges = face.getEdges();
            for (let i = 0; i < edges.length; i++) {
                this.drawEdge(edges[i]);
            }
        }

        // Code for showing boundary cycles
        if (this.showBoundaries) {
            let cycles = this.mesh.getBoundaryCycles();
            for (let i = 0; i < cycles.length; i++) {
                for (let j = 0; j < cycles[i].length; j++) {
                    this.drawEdge(cycles[i][j], [0, 1, 0.5]);
                }
            }
            this.genus = this.mesh.getGenus();
        }
        drawer.repaint(this.camera);
    }
}
