precision mediump float;

#define INF 1.0e+12
#define EPS 1.0e-3 // Reflect/shadow/transmission ray offset
#define MAX_RECURSION 3 // Maximum depth for rays
#define MAX_LIGHTS 10
#define MAX_MATERIALS 10
#define M_PI 3.1415926535897932384626433832795

/*******************************************
                DATA TYPES
********************************************/
struct Material {
  vec3 kd;
  vec3 ks;
  vec3 ka;
  vec3 kt;
  float shininess;
  float refraction;
  int special;
};

struct Light {
    vec3 pos;
    vec3 color;
    vec3 atten;
    vec3 towards;
    float angle;
};

struct Ray {
  vec3 p0;
  vec3 v;
};

struct Intersection {
  vec3 p; // Point of intersection
  vec3 n; // Normal of intersection
  int mIdx; // Index into materials array
  float sCoeff; // Coefficient for checkerboard or special material
};


/*******************************************
                UNIFORMS
********************************************/
// Uniforms set from Javascript that are constant
// over all fragments
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform int numMaterials;
uniform Material materials[MAX_MATERIALS];
uniform int showLights;
uniform float beaconRadius;

// Ray tracer special options
uniform int orthographic;

// Camera parameters
uniform vec3 cameraPos;
uniform vec3 right;
uniform vec3 up;
uniform float fovx;
uniform float fovy; 
// towards is the cross product (up x right)
// right x towrads = up


/*******************************************
           RAY CASTING FUNCTIONS
********************************************/

// TODO: Put helper functions here if you'd like
float triangleArea(vec3 a, vec3 b, vec3 c){
    vec3 v = b - a;
    vec3 u = c - a;
    return length(cross(v,u)) / 2.0;
}

vec3 triangleNormal(vec3 a, vec3 b, vec3 c){
    vec3 v = b - a;
    vec3 u = c - a;
    return cross(v,u);
}

vec3 barycentricCoordinates(vec3 a, vec3 b, vec3 c, vec3 p){
    float fullArea = triangleArea(a,b,c);
    float alpha = triangleArea(b,c,p) / fullArea;
    float beta = triangleArea(a,c,p) / fullArea;
    float gamma = triangleArea(a,b,p) / fullArea;
    return vec3(alpha, beta, gamma);
}

vec2 quadraticFormula(float a, float b, float c){
	vec2 solutions = vec2(INF,INF);

    if((a != 0.0) && ((b*b) - 4.0 * a * c >= 0.0)){
        solutions[0] = ((-1.0 * b) - sqrt((b*b) - 4.0 * a * c)) / (2.0 * a);
	    solutions[1] = ((-1.0 * b) + sqrt((b*b) - 4.0 * a * c)) / (2.0 * a);
    }

	return solutions;
}

float angleRad(vec3 u, vec3 v){
    return acos(dot(normalize(u),normalize(v)));
}
float angleDeg(vec3 u, vec3 v){
    return acos(dot(normalize(u),normalize(v))) * (180.0 / M_PI);
}

/** TODO: PUT YOUR CODE HERE **/


/**
* Given a point on a plane and a normal, intersect a ray
* with the plane they determine
*
* @param {Ray} ray : The ray in world coordinates
* @param {vec3} n : The plane normal
* @param {vec3} p : A point on the plane
* @param {int} mIdx : Array index of material that the plane is made of
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectPlane(Ray ray, vec3 n, vec3 p, int mIdx, out Intersection intersect) {
    float denom = dot(ray.v, n);
    float t = INF;
    if (abs(denom) > 0.0) {
        // The ray is not parallel to the plane
        float num = dot(p - ray.p0, n);
        t = num / denom;
        if (t > 0.0) {
            // Plane is in front of ray
            intersect.p = ray.p0 + t*ray.v;
            intersect.n = n;
            intersect.mIdx = mIdx;
        }
        else {
            t = INF;
        }
    }
    return t;
}


/**
* Intersect a ray with a given triangle /\abc, assuming a, b, and c
* have been specified in CCW order with respect to the triangle normal
*
* @param {Ray} ray : The ray in world coordinates
* @param {vec3} a : Point a on the triangle
* @param {vec3} b : Point b on the triangle
* @param {vec3} c: Point c on the triangle
* @param {int} mIdx : Array index of material that the triangle is made of
* @param {mat4} MInv: Inverse of the transformation M that's applied to the triangle before ray intersection
* @param {mat3} N: The normal transformation associated to M
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectTriangle(Ray ray, vec3 a, vec3 b, vec3 c,
                            int mIdx, mat4 MInv, mat3 N,
                            out Intersection intersect) {
    intersect.mIdx = mIdx; // Store away the material index


/** TODO: PUT YOUR CODE HERE **/
    // TODO: The below three are dummy values
    intersect.p = vec3(0, 0, 0);
    intersect.n = vec3(0, 0, 0);

    Ray OGray = ray;
    vec4 Tp = MInv * vec4(ray.p0,1.0);
    ray.p0 = vec3(Tp[0], Tp[1], Tp[2]);
    vec4 Tv = MInv * vec4(ray.v,0.0);
    ray.v = vec3(Tv[0], Tv[1], Tv[2]);

    float epsilon = 0.001;
    
    intersect.n = normalize(N * triangleNormal(a,b,c));
    float t = rayIntersectPlane(ray, intersect.n, a, mIdx, intersect);
    intersect.p = ray.p0 + (t * ray.v);

    vec3 bCoords = barycentricCoordinates(a,b,c,intersect.p);
    if(bCoords[0] + bCoords[1] + bCoords[2] > 1.0+epsilon || bCoords[0] + bCoords[1] + bCoords[2] < 1.0-epsilon){
        return INF;
    }

    return t;
}


/**
* Intersect a ray with a given sphere
*
* @param {Ray} ray : The ray in world coordinates
* @param {vec3} c : Center of the sphere
* @param {float} r : Radius of the sphere
* @param {int} mIdx : Array index of material that the sphere is made of
* @param {mat4} MInv: Inverse of the transformation M that's applied to the sphere before ray intersection
* @param {mat3} N: The normal transformation associated to M
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectSphere(Ray ray, vec3 c, float r,
                            int mIdx, mat4 MInv, mat3 N,
                            out Intersection intersect) {
    intersect.mIdx = mIdx; // Store away the material index

/** TODO: PUT YOUR CODE HERE **/
    // TODO: The below three are dummy values
    intersect.p = vec3(0, 0, 0);
    intersect.n = vec3(0, 0, 0);
    intersect.sCoeff = 1.0; // TODO: Change this for special material extra task

    // transform rays
    Ray OGray = ray;
    vec4 Tp = MInv * vec4(ray.p0,1.0);
    ray.p0 = vec3(Tp[0], Tp[1], Tp[2]);
    vec4 Tv = MInv * vec4(ray.v,0.0);
    ray.v = vec3(Tv[0], Tv[1], Tv[2]);

    vec3 w = ray.p0 - c;
    float determinant = (dot(ray.v, w) * dot(ray.v, w)) - (dot(ray.v, ray.v) * (dot(w,w) - (r*r))); //(vec3.dot(v, w) ** 2) - (vec3.dot(v,v) * (vec3.dot(w,w) - (r ** 2)));;

    float t1 = (dot(-1.0 * ray.v, w) -  sqrt(determinant)) / dot(ray.v, ray.v); // (vec3.dot(vec3.scale(out2, v, -1), w) - Math.sqrt(determinant)) / vec3.dot(v,v);
    float t2 = (dot(-1.0 * ray.v, w) +  sqrt(determinant)) / dot(ray.v, ray.v);

    
    if(t1 <= t2){
        if(t1 > 0.0){
            intersect.p = ray.p0 + (t1 * ray.v);
            intersect.n = normalize(N * (intersect.p - c));
            intersect.p = OGray.p0 + (t1 * OGray.v);

            return t1;
        }
        else{
            return INF;
        }
    }

    if(t2 < t1){
        if(t2 > 0.0){
            intersect.p = ray.p0 + (t2 * ray.v);
            intersect.n = normalize(N * (intersect.p - c));
            intersect.p = OGray.p0 + (t2 * OGray.v);

            return t2;
        }
        else{
            return INF;
        }
    }

    return INF;
}


/**
* Intersect a ray with a (possibly transformed) box, whose extent
* in untransformed space is [center[0]-width/2, center[0]+width/2],
*                           [center[1]-height/2, center[1]+height/2],
*                           [center[2]-length/2, center[2]+length/2]
*
* @param {Ray} ray : The ray in world coordinates
* @param {float} W : Extent of the box along the x dimension
* @param {float} H : Extent of the box along the y dimension
* @param {float} L : Extent of the box along the z dimension
* @param {vec3} c : Center of the box
* @param {int} mIdx : Array index of material that the box is made of
* @param {mat4} MInv: Inverse of the transformation M that's applied to the box before ray intersection
* @param {mat3} N: The normal transformation associated to M
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectBox(Ray ray, float W, float H, float L,
                        vec3 c, int mIdx, mat4 MInv, mat3 N,
                        out Intersection intersect) {
    intersect.mIdx = mIdx; // Store away the material index

/** TODO: PUT YOUR CODE HERE **/
    // TODO: The below three are dummy values
    intersect.p = vec3(0, 0, 0);
    intersect.n = vec3(0, 0, 0);
    intersect.sCoeff = 1.0; // TODO: Change this for special material extra task

    // transform the box
    Ray OGray = ray;
    vec4 Tp = MInv * vec4(ray.p0,1.0);
    ray.p0 = vec3(Tp[0], Tp[1], Tp[2]);
    vec4 Tv = MInv * vec4(ray.v,0.0);
    ray.v = vec3(Tv[0], Tv[1], Tv[2]);

    float cx = c[0];
    float cy = c[1];
    float cz = c[2];

    // intersections
    Intersection i1;
    Intersection i2;
    Intersection i3;
    Intersection i4;
    Intersection i5;
    Intersection i6;

    // center points of the faces of the box
    vec3 XPoint  = vec3(cx + (W/2.0), cy, cz);
    vec3 NXPoint = vec3(cx - (W/2.0), cy, cz);
    vec3 YPoint  = vec3(cx, cy + (H/2.0), cz);
    vec3 NYPoint = vec3(cx, cy - (H/2.0), cz);
    vec3 ZPoint  = vec3(cx, cy, cz + (L/2.0));
    vec3 NZPoint = vec3(cx, cy, cz - (L/2.0));

    // normals of the faces
    vec3 XNormal  = vec3( 1.0,0.0,0.0); 
    vec3 NXNormal = vec3(-1.0,0.0,0.0);
    vec3 YNormal  = vec3(0.0, 1.0,0.0);
    vec3 NYNormal = vec3(0.0,-1.0,0.0);
    vec3 ZNormal  = vec3(0.0,0.0, 1.0);
    vec3 NZNormal = vec3(0.0,0.0,-1.0);

    // find when the ray hits a plane
    // then check if the intersection point is on a face of the box
    float t = INF;

    float t1 = rayIntersectPlane(ray, XNormal, XPoint, mIdx, i1);
    if((i1.p)[1] < (cy + H/2.0) && (i1.p)[1] > (cy - H/2.0) && (i1.p)[2] < (cz + L/2.0) && (i1.p)[2] > (cz - L/2.0)){
        if (t1 < t){
            t = t1;
            intersect.p = i1.p;
            intersect.n = i1.n;
            intersect.sCoeff = cos(intersect.p[1]*10.0) * cos(intersect.p[2]*10.0);
        }
    }

    float t2 = rayIntersectPlane(ray, NXNormal, NXPoint, mIdx, i2);
    if((i2.p)[1] < (cy + H/2.0) && (i2.p)[1] > (cy - H/2.0) && (i2.p)[2] < (cz + L/2.0) && (i2.p)[2] > (cz - L/2.0)){
        if (t2 < t){
            t = t2;
            intersect.p = i2.p;
            intersect.n = i2.n;
            intersect.sCoeff = cos(intersect.p[1]*10.0) * cos(intersect.p[2]*10.0);
        }
    }

    float t3 = rayIntersectPlane(ray, YNormal, YPoint, mIdx, i3);
    if((i3.p)[0] < (cx + L/2.0) && (i3.p)[0] > (cx - L/2.0) && (i3.p)[2] < (cz + W/2.0) && (i3.p)[2] > (cz - W/2.0)){
        if (t3 < t){
            t = t3;
            intersect.p = i3.p;
            intersect.n = i3.n;
            intersect.sCoeff = cos(intersect.p[0]*10.0) * cos(intersect.p[2]*10.0);
        }
    }

    float t4 = rayIntersectPlane(ray, NYNormal, NYPoint, mIdx, i4);
    if((i4.p)[0] < (cx + L/2.0) && (i4.p)[0] > (cx - L/2.0) && (i4.p)[2] < (cz + W/2.0) && (i4.p)[2] > (cz - W/2.0)){
        if (t4 < t){
            t = t4;
            intersect.p = i4.p;
            intersect.n = i4.n;
            intersect.sCoeff = cos(intersect.p[0]*10.0) * cos(intersect.p[2]*10.0);
        }
    }

    float t5 = rayIntersectPlane(ray, ZNormal, ZPoint, mIdx, i5);
    if((i5.p)[0] < (cx + L/2.0) && (i5.p)[0] > (cx - L/2.0) && (i5.p)[1] < (cy + H/2.0) && (i5.p)[1] > (cy - H/2.0)){
        if (t5 < t){
            t = t5;
            intersect.p = i5.p;
            intersect.n = i5.n;
            intersect.sCoeff = cos(intersect.p[0]*10.0) * cos(intersect.p[1]*10.0);
        }
    }

    float t6 = rayIntersectPlane(ray, NZNormal, NZPoint, mIdx, i6);
    if((i6.p)[0] < (cx + L/2.0) && (i6.p)[0] > (cx - L/2.0) && (i6.p)[1] < (cy + H/2.0) && (i6.p)[1] > (cy - H/2.0)){
        if (t6 < t){
            t = t6;
            intersect.p = i6.p;
            intersect.n = i6.n;
            intersect.sCoeff = cos(intersect.p[0]*10.0) * cos(intersect.p[1]*10.0);
        }
    }

    if(intersect.sCoeff >= 0.0){
        intersect.sCoeff = 1.0;
    }
    else{
        intersect.sCoeff = 0.0;
    }
    intersect.p = OGray.p0 + (t * OGray.v);
    intersect.n = normalize(N * intersect.n);
    
    return t;
    }


/**
* Intersect a ray with a given cylinder
*
* @param {Ray} ray : The ray in world coordinates
* @param {vec3} c : Center of cylinder
* @param {float} r : Radius of cylinder
* @param {float} h : Height of cylinder
* @param {int} mIdx : Array index of material that the cylinder is made of
* @param {mat4} MInv: Inverse of the transformation M that's applied to the cylinder before ray intersection
* @param {mat3} N: The normal transformation associated to M
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectCylinder(Ray ray, vec3 c, float r, float h,
                            int mIdx, mat4 MInv, mat3 N,
                            out Intersection intersect) {
    intersect.mIdx = mIdx; // Store away the material index
    intersect.sCoeff = 1.0; // TODO: Change this for special material extra task
/** TODO: PUT YOUR CODE HERE **/
    // TODO: The below three are dummy values

    Ray OGray = ray;
    vec4 Tp = MInv * vec4(ray.p0,1.0);
    ray.p0 = vec3(Tp[0], Tp[1], Tp[2]);
    vec4 Tv = MInv * vec4(ray.v,0.0);
    ray.v = vec3(Tv[0], Tv[1], Tv[2]);
    
    float t = INF;
    float cy = c[1];
    float rp0y = (ray.p0)[1];
    vec3 topCenter = vec3(c[0], c[1] + h/2.0, c[2]);
    vec3 bottomCenter = vec3(c[0], c[1] - h/2.0, c[2]);
    vec3 up = vec3(0,1,0);
    vec3 down = vec3(0,-1,0);
    Intersection i;

    //float rayIntersectPlane(Ray ray, vec3 n, vec3 p, int mIdx, out Intersection intersect) 
    if(rp0y > cy + h/2.0){
        // check if ray hit the bottom
        t = rayIntersectPlane(ray, up, topCenter, mIdx, i); 
        if(distance(topCenter, ray.p0 + t * ray.v) <= r){
            intersect.p = ray.p0 + t * ray.v;
            intersect.n = up;
            return t;
        }
    }
    else if(rp0y < cy - h/2.0){
        // check if ray hit the top
        t = rayIntersectPlane(ray, down, bottomCenter, mIdx, i); 
        if(distance(bottomCenter, ray.p0 + t * ray.v) <= r){
            intersect.p = ray.p0 + t * ray.v;
            intersect.n = down;
            return t;
        }
    }

    vec3 cFlat = vec3(c[0], (ray.p0)[1], c[2]);
    vec3 vFlat = vec3((ray.v)[0], 0.0, (ray.v)[2]);
    vec3 w = ray.p0 - cFlat;

    float a = dot(vFlat,vFlat);
    float b = dot(2.0 * w, vFlat);
    float cc = dot(w,w) - (r*r); 

    vec2 times = quadraticFormula(a,b,cc);
    t = min(times[0], times[1]);
    if(t < 0.0){
        t = INF;
    }

    intersect.p = OGray.p0 + (t * OGray.v); //ray.p0 + (t * ray.v);
    intersect.n = normalize(N * ((ray.p0 + (t * vFlat)) - cFlat));

    if((intersect.p)[1] > topCenter[1] || (intersect.p)[1] < bottomCenter[1]){
        t = INF;
    }

    return t;
}


/**
* Intersect a ray with a given cone
*
* @param {Ray} ray : The ray in world coordinates
* @param {vec3} c : Center of cone
* @param {float} r : Radius of cone
* @param {float} h : Height of cone
* @param {int} mIdx : Array index of material that the cone is made of
* @param {mat4} MInv: Inverse of the transformation M that's applied to the cone before ray intersection
* @param {mat3} N: The normal transformation associated to M
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectCone(Ray ray, vec3 c, float r, float h,
                            int mIdx, mat4 MInv, mat3 N,
                            out Intersection intersect) {
    intersect.mIdx = mIdx; // Store away the material index
    intersect.sCoeff = 1.0; // TODO: Change this for special material extra task
/** TODO: PUT YOUR CODE HERE **/
    // TODO: The below three are dummy values
/**
    Ray OGray = ray;
    vec4 Tp = MInv * vec4(ray.p0,1.0);
    ray.p0 = vec3(Tp[0], Tp[1], Tp[2]);
    vec4 Tv = MInv * vec4(ray.v,0.0);
    ray.v = vec3(Tv[0], Tv[1], Tv[2]);

    float t = INF;

    float cx = c[0];
    float cy = c[1];
    float cz = c[2];

    float p0x = ray.p0[0];
    float p0y = ray.p0[1];
    float p0z = ray.p0[2];

    float vx = ray.v[0];
    float vy = ray.v[1];
    float vz = ray.v[2];

    float k = r/h;
    float yk = p0y - (cy + h);

    vec3 down = vec3(0.0,-1.0,0.0);

    // draw circle
    if(p0y < cy){
        t = rayIntersectPlane(ray, down, c, mIdx, intersect); 
        if(distance(c, ray.p0 + (t * ray.v)) <= r){
            intersect.p = ray.p0 + (t * ray.v);
            intersect.n = down;
            return t;
        }
    }

    float a = vx*vx + vz*vz - k*k*vy*vy;
    float b = 2.0 * (vx*p0x + vz*p0z - k*k*yk*vy);
    float cc = p0x*p0x + p0z*p0z - k*k*yk*yk;

    vec2 times = quadraticFormula(a,b,cc);
    t = min(times[0], times[1]);
    if(t < 0.0){
        t = INF;
    }
    
    intersect.p = OGray.p0 + (t * OGray.v);
    vec3 nFlat = intersect.p - c;
    intersect.n = normalize( N * vec3(nFlat[0], length(vec3(nFlat[0], 0.0, nFlat[2])) * tan(k * (M_PI / 180.0)), nFlat[2]));	

    if((intersect.p)[1] < cy || (intersect.p)[1] > cy + h){
        t = INF;
    }

    return t;
    */
    return INF;
}


/**
* A function which intersects a ray with a scene, returning the
* t parameter of the closest intersection, or INF if no intersection
* happened, along with an out parameter storing the point, normal,
* and material of the intersection
* NOTE: This function is merely declared here; it is defined in its
* entirety in Javascript before this shader is compiled, since information
* about the scene must be hardcoded into the shader
*
* @param {Ray} ray : The ray in world coordinates
* @param {Intersection (out)} intersect : The intersection
*
* @returns {float} t : Parameter t so that point of intersection is ray.p0 + t*ray.v
*/
float rayIntersectScene(Ray ray, out Intersection intersect){return INF;}


/*******************************************
        RAY ILLUMINATION FUNCTIONS
********************************************/

/**
* Pull a material out of the list of materials, based on its
* index.  This function is necessary, because it is not possible
* to use non-constant indices into arrays in GLSL, so one must
* loop over the entire array of materials to find the right one
*
* @param {int} mIdx : The index into the materials array of the
*                     material we seekd
*
* @returns {Material} m : The appropriate material struct
*/
Material getMaterial(int mIdx) {
    Material m;
    for (int i = 0; i < MAX_MATERIALS; i++) {
        if (i == mIdx) {
            m = materials[i];
        }
        if (i >= numMaterials) {
            break;
        }
    }
    return m;
}

/**
* Determine whether a point is in the shadow of a light
*
* @param {Intersection} intersect : Intersection point we're checking
* @param {int} lightIndex : Index into the array of lights of
*                           the light we want to check
*/
bool pointInShadow(Intersection intersect, Light l) {

/** TODO: PUT YOUR CODE HERE **/
    float epsilon = 0.001;
    Ray reverseRay;
    reverseRay.p0 = intersect.p + (epsilon * (l.pos - intersect.p));
    reverseRay.v = l.pos - reverseRay.p0;
    //float rayIntersectScene(Ray ray, out Intersection intersect){return INF;}
    float t = rayIntersectScene(reverseRay, intersect);
    if(t < 1.0 && t > 0.0){
        return true;
    }
    return false;
}

/**
* Get the phong illumination color
*/
vec3 getPhongColor(Intersection intersect, Material m, vec3 eye) {
    vec3 zeroVec = vec3(0.0, 0.0, 0.0);
    vec3 color = zeroVec;
    vec3 colour = zeroVec;
    vec3 ci = zeroVec;
    vec3 unitLight = zeroVec;
    vec3 diffuse = zeroVec;
    vec3 specular = zeroVec;
    vec3 lightVec = zeroVec;
    float d = 0.0;
    float specCoeff = 0.0;
    float diffCoeff = 0.0;
    bool shadowed;
    Light lamp;
    // To help with debugging, color the fragment based on the
    // normal of the intersection.  But this should eventually
    // be replaced with code to do Phong illumination below
    //color = 0.5*(intersect.n + 1.0);

/** TODO: PUT YOUR CODE HERE **/
    for(int l = 0; l < MAX_LIGHTS; l++){
        if(l < numLights){
            
            lamp = lights[l];
            shadowed = pointInShadow(intersect, lamp);

            lightVec = lamp.pos - intersect.p;
            d = length(lightVec);
            unitLight = normalize(lightVec);

            diffCoeff = dot(intersect.n, unitLight);
            specCoeff = dot(normalize(intersect.p - eye), reflect(unitLight, intersect.n));
            
            if(diffCoeff < 0.0){
                diffCoeff = 0.0;
            }
            if(specCoeff < 0.0){
                specCoeff = 0.0;
            }

            if(m.special == 1){
                m.kd *= intersect.sCoeff;
            }
            diffuse = m.kd * diffCoeff;
            
            specular = m.ks * pow(specCoeff,m.shininess);

            ci = lamp.color / (lamp.atten[0] + lamp.atten[1] * d + lamp.atten[2] * (d * d));
            colour = ci * (diffuse + specular);

            if(shadowed){
                colour = zeroVec;
            }

            if(angleRad(lamp.towards, -1.0 * lightVec) > lamp.angle){
                colour = zeroVec;
            }

            color += colour;
        }
        else{
            break;
        }
    }

    return color;
}


/**
*
*/
varying vec2 v_position;
Ray getRay() {
    Ray ray;
    ray.p0 = cameraPos;
    // TODO: Finish constructing ray by figuring out direction, using
    // v_position.x, v_position.y, fovx, fovy, up, and right
    if (orthographic == 1) {
        // TODO: Fill in code for constructing orthographic rays
        // (You can ignore this if you aren't doing the orthographic extra task)

        /** TODO: PUT YOUR CODE HERE **/
        vec3 towards = cross(up,right);
        ray.p0 = cameraPos + right * v_position.x + up * v_position.y;
        ray.v = towards;
    }
    else {
        // TODO: Fill in ordinary perspective ray based on fovx and fovy (the default option)
        
        /** TODO: PUT YOUR CODE HERE **/
        vec3 towards = cross(up,right);
        ray.v = towards + (right * v_position.x * tan(fovx/2.0)) + (up * v_position.y * tan(fovy/2.0));
    }
    return ray;
}

void showLightBeacons(Ray rayInitial, float tInitial) {
    // Show light beacons if the user so chooses
    // (NOTE: This requires a working implementation of rayIntersectSphere)
    mat4 identity4 = mat4(1.0);
    mat3 identity3 = mat3(1.0);
    Intersection intersect;
    if (showLights == 1) {
        for (int i = 0; i < MAX_LIGHTS; i++) {
            if (i < numLights) {
                Light light = lights[i];
                float tlight = rayIntersectSphere(rayInitial, light.pos, beaconRadius,
                                                  0, identity4, identity3, intersect);
                if (tlight < tInitial) {
                    gl_FragColor = vec4(light.color, 1.0);
                }
            }
        }
    }
}

void main() {
    Ray ray = getRay();
    Ray rayInitial = ray;
    bool insideObj = false;
    Intersection intersect;
    intersect.sCoeff = 1.0;
    vec3 color = vec3(0.0, 0.0, 0.0);
    vec3 weight = vec3(1.0, 1.0, 1.0);
    float t;
    float tInitial;
    float epsilon = 0.001;
    for (int depth = 0; depth < MAX_RECURSION; depth++) {
        t = rayIntersectScene(ray, intersect);
        if (depth == 0) {
            tInitial = t;
        }
        if (t < INF) {
            Material m = getMaterial(intersect.mIdx);
            // Figure out whether the ray is inside the object it
            // intersected by using the dot product between a vector
            // from the endpoint of the ray and the intersection
            // point and the intersection normal
            if (dot(ray.p0 - intersect.p, intersect.n) < 0.0) {
                intersect.n *= -1.0;
                insideObj = true;
            }
            else {
                insideObj = false;
            }
            color += weight*getPhongColor(intersect, m, rayInitial.p0);
            // TODO: Reflect ray, multiply weight by specular of this object,
            // and recursively continue
            // If doing extra task on transmission, only reflect if the
            // transmission coefficient kt is zero in all components
            // Otherwise, do transmission with snell's law
            if(m.ks != vec3(0.0, 0.0, 0.0)){
                ray.v = reflect(ray.v, intersect.n);
                ray.p0 = ray.p0 + (epsilon * ray.v);
                weight *= m.ks;
            }

/** TODO: PUT YOUR CODE HERE **/
        }
        else {
            // Ray doesn't intersect anything, so no use continuing
            break;
        }
    }
    gl_FragColor = vec4(color, 1.0);
    showLightBeacons(rayInitial, tInitial);
}
