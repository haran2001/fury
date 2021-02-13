/* SDF fragment shader declaration */

//VTK::ValuePass::Dec
in vec4 vertexMCVSOutput;

in vec3 centerWCVSOutput;
flat in int primitiveVSOutput;
in float scaleVSOutput;
in vec3 directionVSOutput;

uniform mat4 MCVCMatrix;
uniform mat4 MCWCMatrix;
uniform mat4 WCVCMatrix;
uniform mat4 WCMCMatrix;
uniform mat4 VCMCMatrix;
uniform float u_time;

#define PI 3.14159265

struct Superellipsoid
{
	vec3 Center;
	vec3 Radius;
    vec2 Exponent;
    mat3 Orientation;
};

mat4 QuaternionToMatrix(vec4 q)
{
    mat4 m1 = mat4( q.w,  q.z, -q.y, q.x,
                   -q.z,  q.w,  q.x, q.y,
                    q.y, -q.x,  q.w, q.z,
                   -q.x, -q.y, -q.z, q.w);

    mat4 m2 = mat4( q.w,  q.z, -q.y, -q.x,
                   -q.z,  q.w,  q.x, -q.y,
                    q.y, -q.x,  q.w, -q.z,
                    q.x,  q.y,  q.z,  q.w);

    mat4 m = m1 * m2;

    return(m);
}

mat4 AxisAngleToMatrix(vec3 axis, float angle)
{
    float s = sin(angle / 2.0);
    float c = cos(angle / 2.0);

    vec4 q = vec4(s, s, s, c);
    q.x *= axis.x;
    q.y *= axis.y;
    q.z *= axis.z;

    mat4 m = QuaternionToMatrix(q);
    return(m);
}


mat4 rotationAxisAngle( vec3 v, float angle )
{
    float s = sin(angle);
    float c = cos(angle);
    float ic = 1.0 - c;

    return mat4( v.x*v.x*ic + c,     v.y*v.x*ic - s*v.z, v.z*v.x*ic + s*v.y, 0.0,
                 v.x*v.y*ic + s*v.z, v.y*v.y*ic + c,     v.z*v.y*ic - s*v.x, 0.0,
                 v.x*v.z*ic - s*v.y, v.y*v.z*ic + s*v.x, v.z*v.z*ic + c,     0.0,
                 0.0,                0.0,                0.0,                1.0 );
}


mat4 translate( float x, float y, float z )
{
    return mat4( 1.0, 0.0, 0.0, 0.0,
                 0.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 0.0,
                 x,   y,   z,   1.0 );
}


float sdEllipsoid( vec3 p, vec3 r )
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
}

float sdTest( vec3 pos , Superellipsoid se )
{
    vec3 e = vec3(vec2(1.0) / se.Exponent.xy, se.Exponent.x / se.Exponent.y);
    vec3 g = 2.0 * e; // prep exponents
    vec3 invr = vec3(1.0) / se.Radius;
    vec3 p = pos + 1.0 * vec3(0.0);
    vec3 A = p * invr; // coeff for x and y
    vec3 B = pow(A * A, e.xxy); //prep x and y with exponents
    float E = B.x + B.y;
    float F = pow(E, e.z);
    float P = F + B.z;

    float K = pow(P, se.Exponent.y) - 1.0;
    return(K);

}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}


float sdTorus( vec3 p, vec2 t )
{
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}

float sdCapsule(vec3 p, vec3 a, vec3 b, float r){
    vec3 ab = vec3(b - a);
    vec3 ap = vec3(p - a);

    float t = clamp(dot(ab, ap) / dot(ab, ab), 0.0, 1.0);
    vec3 c = a + t*ab;
    return length(p-c) - r;
}

float map( in vec3 position )
{

    mat4 rot = rotationAxisAngle( normalize(directionVSOutput), 90.0 );
    mat4 tra = translate( 0.0, 1.0, 0.0 );
    mat4 txi = tra * rot; 

    vec3 pos = (txi*vec4(position  - centerWCVSOutput, 0.0)).xyz;
	
    float d1;
	
    if(primitiveVSOutput==1){
		// d1 = sdSphere((pos)/scaleVSOutput, 0.25)*scaleVSOutput;
		d1 = sdSphere((pos)/scaleVSOutput, 0.25)*scaleVSOutput;
    }
    
    else if(primitiveVSOutput==2){
    	d1 = sdTorus((pos)/scaleVSOutput, vec2(0.4, 0.1))*scaleVSOutput;
    }
    
    else if(primitiveVSOutput==3){
        d1 = sdEllipsoid((pos)/scaleVSOutput, vec3(0.1, 0.1, 0.3))*scaleVSOutput;
    }

     else if(primitiveVSOutput==4){
        vec3 center = vec3(0.0, 0.0, 0.0);
        // vec3 radius = vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0);
        // vec3 radius = vec3(1.0 / 12.0, 1.0 / 12.0, 1.0 / 8.0);
        // vec3 radius = vec3(1.0, 1.0, 1.5);
        float a = 0.5;
        float b = 0.5;
        vec3 radius = vec3(a , a, b);
        vec2 e = vec2(1.0);
        
        vec3 axis1 = vec3(0.0, 1.0, 1.0);
        vec3 axis = axis1;
        mat4 o = AxisAngleToMatrix(normalize(axis), PI);

        // Superellipsoid se = {center, radius, e, mat3(o)};
        Superellipsoid se;
        se.Center = center;
        se.Radius = radius;
        se.Exponent = e;
        se.Orientation = mat3(o);
        
        // mat3 invm = transpose(se.Orientation);
        // pos = invm * (pos- se.Center);
        vec3 shift = vec3(abs(sin(u_time)));
        d1 = sdTest((pos - shift)/scaleVSOutput, se)*scaleVSOutput;
    }

     else if(primitiveVSOutput==5){
        //  vec3 shift1 = vec3(0.0, -0.2, -0.2);
         vec3 shift1 = vec3(0.0, -0.1, -0.1);
         vec3 shift2 = vec3(0.0);
         vec3 a = vec3(0.0, 0.1, 0.6);
         vec3 b = vec3(0.0, 0.2, 0.6);
         float t = 0.2;
        d1 = sdCapsule((pos - shift1)/scaleVSOutput, a - shift2, b - shift2, t)*scaleVSOutput;
    }

    return d1;
}


vec3 calculateNormal( in vec3 position )
{
    vec2 e = vec2(0.001, 0.0);
    return normalize( vec3( map(position + e.xyy) - map(position - e.xyy),
    						map(position + e.yxy) - map(position - e.yxy),
    						map(position + e.yyx) - map(position - e.yyx)
                          )
                    );

}


float castRay( in vec3 ro, vec3 rd )
{
    float t = 0.0;
    for(int i=0; i < 4000; i++){

    	vec3 position = ro + t * rd;
    	float  h = map(position);
    	t += h;

    	if ( t > 20.0 || h < 0.001) break;
    }
    return t;
}

