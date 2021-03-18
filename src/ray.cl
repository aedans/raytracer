/**
@file

Implements Mersenne twister generator. 

M. Matsumoto, T. Nishimura, Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator, ACM Transactions on Modeling and Computer Simulation (TOMACS) 8 (1) (1998) 3–30.
*/
#define MT19937_FLOAT_MULTI 2.3283064365386962890625e-10f
#define MT19937_DOUBLE2_MULTI 2.3283064365386962890625e-10
#define MT19937_DOUBLE_MULTI 5.4210108624275221700372640e-20

#define MT19937_N 624
#define MT19937_M 397
#define MT19937_MATRIX_A 0x9908b0df   /* constant vector a */
#define MT19937_UPPER_MASK 0x80000000 /* most significant w-r bits */
#define MT19937_LOWER_MASK 0x7fffffff /* least significant r bits */

/**
State of MT19937 RNG.
*/
typedef struct{
	uint mt[MT19937_N]; /* the array for the state vector  */
	int mti;
} mt19937_state;

/**
Generates a random 32-bit unsigned integer using MT19937 RNG.

@param state State of the RNG to use.
*/
#define mt19937_uint(state) _mt19937_uint(&state)
uint _mt19937_uint(mt19937_state* state){
    uint y;
    uint mag01[2]={0x0, MT19937_MATRIX_A};
    /* mag01[x] = x * MT19937_MATRIX_A  for x=0,1 */
	
	if(state->mti<MT19937_N-MT19937_M){
		y = (state->mt[state->mti]&MT19937_UPPER_MASK)|(state->mt[state->mti+1]&MT19937_LOWER_MASK);
		state->mt[state->mti] = state->mt[state->mti+MT19937_M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	else if(state->mti<MT19937_N-1){
		y = (state->mt[state->mti]&MT19937_UPPER_MASK)|(state->mt[state->mti+1]&MT19937_LOWER_MASK);
		state->mt[state->mti] = state->mt[state->mti+(MT19937_M-MT19937_N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	else{
        y = (state->mt[MT19937_N-1]&MT19937_UPPER_MASK)|(state->mt[0]&MT19937_LOWER_MASK);
        state->mt[MT19937_N-1] = state->mt[MT19937_M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        state->mti = 0;
	}
    y = state->mt[state->mti++];
		
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return y;
}
/**
Seeds MT19937 RNG.

@param state Variable, that holds state of the generator to be seeded.
@param seed Value used for seeding. Should be randomly generated for each instance of generator (thread).
*/
void mt19937_seed(mt19937_state* state, uint s){
    state->mt[0]= s;
	uint mti;
    for (mti=1; mti<MT19937_N; mti++) {
        state->mt[mti] = 1812433253 * (state->mt[mti-1] ^ (state->mt[mti-1] >> 30)) + mti;
		
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt19937[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
    }
	state->mti=mti;
}

/**
Generates a random float using MT19937 RNG.

@param state State of the RNG to use.
*/
#define rand(state) (_mt19937_uint(state)*MT19937_FLOAT_MULTI)

/* Code starts here */

#define EPSILON 1e-4

struct __attribute__((packed)) ray {
  float3 origin;
  float3 direction;
};

struct __attribute__((packed)) material {
  float mirror_coeff;
  float transparent_coeff;
  float3 diffuse_color;  
};

struct __attribute__((packed)) light {
  float3 color;
  float3 point;
  float radius;
};

struct __attribute__((packed)) triangle {
  float3 a;
  float3 b;
  float3 c;
  float3 anormal;
  float3 bnormal;
  float3 cnormal;
  float2 auv;
  float2 buv;
  float2 cuv;
  int material;
};

struct __attribute__((packed)) bounds {
  float3 min;
  float3 max;
  int leaf;
  int a;
  int b;
};

struct hitrec {
  float t;
  float3 intersection;
  float3 normal;
  float2 uv;
  int material;
};

bool intersect_triangle(
  struct triangle triangle,
  struct ray ray,
  float tmin,
  float tmax,
  struct hitrec *hitrec
) {
  float3 p = ray.origin;
  float3 d = normalize(ray.direction);

  float3 a = triangle.a;
  float3 b = triangle.b;
  float3 c = triangle.c;

  float3 e1 = b - a;
  float3 e2 = c - a;

  float3 q = cross(d, e2);

  float detA = dot(q, e1);

  if (-EPSILON < detA && detA < EPSILON) {
    return false;
  }

  float3 s = p - a;
  float detB = dot(q, s);
  float B = detB / detA;

  if (!(-EPSILON < B && B < 1 + EPSILON)) {
    return false;
  }

  float3 r = cross(s, e1);
  float detC = dot(r, d);
  float C = detC / detA;

  if (!(-EPSILON < C && C < 1 + EPSILON)) {
    return false;
  }

  float A = 1 - B - C;

  if (!(-EPSILON < A && A < 1 + EPSILON)) {
    return false;
  }

  float detT = dot(r, e2);
  float t = detT / detA;

  if (!(tmin < t && t < tmax)) {
    return false;
  }

  float3 intersection = p + t * d;

  float3 normal = A * triangle.anormal + B * triangle.bnormal + C * triangle.cnormal;
  float2 uv = A * triangle.auv + B * triangle.buv + C * triangle.cuv;

  hitrec->t = t;
  hitrec->intersection = intersection;
  hitrec->normal = normal;
  hitrec->uv = uv;
  hitrec->material = triangle.material;

  return true;
}

bool intersect_triangles(
  __global struct triangle *triangles,
  int num_triangles,
  struct ray ray,
  float tmin,
  float tmax,
  struct hitrec *hitrec
) {
  float min = INFINITY;
  struct hitrec h;

  for (int i = 0; i < num_triangles; ++i) {
    if (intersect_triangle(triangles[i], ray, tmin, tmax, &h) && h.t < min) {
      hitrec->t = h.t;
      hitrec->intersection = h.intersection;
      hitrec->normal = h.normal;
      hitrec->uv = h.uv;
      hitrec->material = h.material;

      min = h.t;
    }
  }

  return min != INFINITY;
}

bool intersect_box(float3 box_min, float3 box_max, struct ray ray) {
  float3 invdir = 1 / ray.direction;

  float tx1 = (box_min.x - ray.origin.x) * invdir.x;
  float tx2 = (box_max.x - ray.origin.x) * invdir.x;

  float tmin = min(tx1, tx2);
  float tmax = max(tx1, tx2);

  float ty1 = (box_min.y - ray.origin.y) * invdir.y;
  float ty2 = (box_max.y - ray.origin.y) * invdir.y;

  tmin = max(tmin, min(ty1, ty2));
  tmax = min(tmax, max(ty1, ty2));

  float tz1 = (box_min.z - ray.origin.z) * invdir.z;
  float tz2 = (box_max.z - ray.origin.z) * invdir.z;

  tmin = max(tmin, min(tz1, tz2));
  tmax = min(tmax, max(tz1, tz2));

  return tmax >= tmin;
}

bool intersect_bounds(
  __global struct triangle *triangles,
  int num_triangles,
  __global struct bounds *bounds,
  int num_bounds,
  struct bounds bound,
  struct ray ray,
  float tmin,
  float tmax,
  struct hitrec *hitrec
) {
  int j = 1;
  struct bounds test[256];
  test[0] = bound;
  
  float min = INFINITY;
  struct hitrec h;

  for (int i = 0; i < j; i++) {
    if (!intersect_box(test[i].min, test[i].max, ray)) {
      continue;
    }

    if (!test[i].leaf) {
      test[j++] = bounds[test[i].a];
      test[j++] = bounds[test[i].b];
    } else {
      if (intersect_triangle(triangles[test[i].a], ray, tmin, tmax, &h) && h.t < min) {
        hitrec->t = h.t;
        hitrec->intersection = h.intersection;
        hitrec->normal = h.normal;
        hitrec->uv = h.uv;
        hitrec->material = h.material;

        min = h.t;
      }

      if (intersect_triangle(triangles[test[i].b], ray, tmin, tmax, &h) && h.t < min) {
        hitrec->t = h.t;
        hitrec->intersection = h.intersection;
        hitrec->normal = h.normal;
        hitrec->uv = h.uv;
        hitrec->material = h.material;

        min = h.t;
      }
    }
  }

  return min != INFINITY;
}

float3 get_diffuse(__global struct material *materials, struct hitrec hitrec){
  return materials[hitrec.material].diffuse_color;
}

float shadow_percent(
  __global struct material *materials,
  __global struct triangle *triangles, 
  int num_triangles, 
  __global struct bounds *bounds,
  int num_bounds,
  float3 lp, 
  float3 point
){
  float dist = distance(lp, point);

  struct ray ray = {point, normalize(lp - point)};
  struct hitrec hitrec;
  bool intersects = intersect_bounds(triangles, num_triangles, bounds, num_bounds, bounds[num_bounds - 1], ray, EPSILON, dist, &hitrec);

  if(!intersects){
    return 1.0;
  }else if(materials[hitrec.material].transparent_coeff < EPSILON){
    return 0.0;
  }
}

float3 traceray(
  __global struct material *materials,
  __global struct light *lights,
  int num_lights,
  __global struct triangle *triangles,
  int num_triangles,
  __global struct bounds *bounds,
  int num_bounds,
  float3 background,
  struct ray ray,
  float tmin,
  float tmax,
  int rec_depth,
  mt19937_state *seed
) {
  struct hitrec hitrec;
  if (!intersect_bounds(triangles, num_triangles, bounds, num_bounds, bounds[num_bounds - 1], ray, tmin, tmax, &hitrec)) {
    return background;
  }

  float3 sample = normalize((float3)(rand(seed), rand(seed), rand(seed)) - 0.5f);
  float d = dot(sample, hitrec.normal);
  if (d < 0) {
    sample = -sample;
    d = -d;
  }

  float3 kd = get_diffuse(materials, hitrec);
  float3 color = (float3)(0.0, 0.0, 0.0);

  for (int i = 0; i < num_lights; ++i) {
    struct light light = lights[i];
    float3 lp = light.point + ((float3)(rand(seed), rand(seed), rand(seed)) - 0.5f) * light.radius;
    float3 c = kd * light.color;

    color += shadow_percent(materials, triangles, num_triangles, bounds, num_bounds, lp, hitrec.intersection) * d * c * max(0.0f, dot(hitrec.normal, normalize(lp - hitrec.intersection)));
  }

  if(rec_depth <= 0){
    return color;
  }

  float mirror_coeff = materials[hitrec.material].mirror_coeff;
  if(mirror_coeff > 0){
    float3 r = ray.direction + 2 * dot(-ray.direction, hitrec.normal) * hitrec.normal;
    ray.origin = hitrec.intersection;
    ray.direction = r;
    float3 reflected_color = traceray(materials, lights, num_lights, triangles, num_triangles, bounds, num_bounds, background, ray, EPSILON, INFINITY, rec_depth - 1, seed);
    color = mirror_coeff * reflected_color + (1 - mirror_coeff) * color;
  }
  
  struct ray sample_ray = {.origin = hitrec.intersection, .direction = sample};
  float3 ray_color = traceray(materials, lights, num_lights, triangles, num_triangles, bounds, num_bounds, background, sample_ray, EPSILON, INFINITY, 0, seed);
  color = (color + d * kd * ray_color) / 2;

  return color;
}

__kernel void ray_kernel(
  __global uint *seeds,
  __global struct ray *rays,
  __global struct material *materials,
  __global struct light *lights,
  int num_lights,
  __global struct triangle *triangles,
  int num_triangles,
  __global struct bounds *bounds,
  int num_bounds,
  float3 background,
  __global float3 *canvas
) {
  int gid = get_global_id(0);
  mt19937_state seed;
  mt19937_seed(&seed, seeds[gid]);
  canvas[gid] = traceray(materials, lights, num_lights, triangles, num_triangles, bounds, num_bounds, background, rays[gid], 1, INFINITY, 1, &seed);
}
