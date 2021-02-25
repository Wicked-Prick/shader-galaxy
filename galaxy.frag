#define PI 3.14159265358979323846 
#define HASHSCALE1 .1031

/* Color palette */
#define iterations 13
#define formuparam 0.53
#define volsteps 13
#define stepsize 0.1
#define zoom   0.800
#define tile   0.850
#define speed  0.010 
#define brightness 0.0035
#define darkmatter 0.300
#define distfading 0.800
#define saturation .990

//uniform sampler2D space; // /home/sabaziy/Документы/glsl/0.jpg
//uniform vec2 iResolution; // size of the preview

//uniform float u_time = iTime;; // clock in seconds
//#iChannel0 'file://home/sabazius/Документы/shaders/noise.jpg'
float Speed = .5;
vec4 MilkyWay;
vec4 BlackHole;

float hash11(float p)
{
	  vec3 p3  = fract(vec3(p) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}
float smoothNoise13( in vec3 x )
{
	vec3 p  = floor(x);
	vec3 f  = smoothstep(0.0, 1.0, fract(x));
	float n = p.x + p.y*57.0 + 113.0*p.z;

	return	mix(
        		mix(
                    mix( hash11( n + 0.0 ), hash11( n + 1.0 ), f.x ),
					mix( hash11( n + 57.0 ), hash11( n + 58.0 ), f.x ),
                    f.y ),
				mix(
                    mix( hash11( n + 113.0 ), hash11( n + 114.0 ), f.x ),
					mix( hash11( n + 170.0 ), hash11( n + 171.0 ), f.x),
                    f.y ),
        		f.z );
}

mat3 ma = mat3( 0.00,  1.60,  1.20, -1.60,  0.72, -0.96, -1.20, -0.96,  1.28 );

float FractionalBrownianMotion( vec3 p )
{
	float f = 0.5000 * smoothNoise13( p );
    p = ma * p * 1.2;
	f += 0.2500 * smoothNoise13( p );
    p = ma * p * 1.3;
	f += 0.1666 * smoothNoise13( p );
    p = ma * p * 1.4;
	f += 0.0834 * smoothNoise13( p );
	return f;
}


float hash( float n )
{
  return fract(cos(n)*41415.92653);
}
mat2 m = mat2( 0.00,  1.60,  1.20, -1.60);

// 2d noise function
float noise( in vec2 x )
{
  vec2 p  = floor(x);
  vec2 f  = smoothstep(0.0, 1.0, fract(x));
  float n = p.x + p.y*57.0;

  return mix(mix( hash(n), hash(n+  1.0),f.x),
    mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y);
}

float fbm( vec2 p )
{
  float f = 0., a = 1., s=0.;
  f += a*noise( p ); p = m*p*1.149; s += a; a *= .75;
  f += a*noise( p ); p = m*p*1.41; s += a; a *= .75;
  f += a*noise( p ); p = m*p*1.51; s += a; a *= .65;
  f += a*noise( p ); p = m*p*1.21; s += a; a *= .35;
  return f/s;
}

vec3 kali_stars(vec3 from, vec3 dir)
{
	//volumetric rendering
	float s=0.1,fade=1.;
	vec3 v=vec3(0.);
	for (int r=0; r<volsteps; r++) {
		vec3 p=from+s*dir*.5;
		p = abs(vec3(tile)-mod(p,vec3(tile*2.))); // tiling fold
		float pa,a=pa=0.;
		for (int i=0; i<iterations; i++) { 
			p=abs(p)/dot(p,p)-formuparam; // the magic formula
			a+=abs(length(p)-pa); // absolute sum of average change
			pa=length(p);
		}
		float dm=max(0.,darkmatter-a*a*.001); //dark matter
		a*=a*a; // add contrast
		if (r>6) fade*=1.-dm; // dark matter, don't render near
		//v+=vec3(dm,dm*.5,0.);
		v+=fade;
		v+=vec3(s,s*s,s*s*s*s)*a*brightness*fade; // coloring based on distance
		fade*=distfading; // distance fading
		s+=stepsize;
	}
    return mix(vec3(length(v)),v,saturation)*.01;
}

vec4 make_MilkyWay(vec2 uv, float rho)
{
  vec4 color;
  float cost = cos(iTime *Speed - 2.*PI * log(rho));
  float sint = sin(iTime *Speed - 2.*PI * log(rho));
  mat2 ma = mat2(cost, -sint, sint, cost);
  uv = uv * ma;
  
  float phi = atan(uv.x, uv.y);
  float co = cos(phi);
  float si = sin(phi);

  vec3 vFbmInput = vec3( uv.x + iTime*0.08,uv.y , 0.0 );
  vec3 vFogColor = vec3(0.78, 0.71, 0.96); 
  float cloud;
  cloud = FractionalBrownianMotion( vFbmInput );
  
  float c = smoothstep(cloud*co + .55 , cloud*si + 1.8 , 1.);

  color = vec4(mix(c * vFogColor ,vec3(0.), vec3(0.)), 1.);
  
  return color;
}

vec4 MakeBlackHole(vec2 uv, float rho)
{
  float waveStrength = .02;
  float frequency = 30.0;
  float waveSpeed = 4.0;
  vec4 sunlightColor = vec4(1.0,0.91,0.75, 1.0);
  float sunlightStrength = 6.0;
  float centerLight = 2.;
  //float oblique = 1.25; 
      
  vec2 center = vec2(0., 0.);

  float t = iTime * waveSpeed;
  float aspectRatio = iResolution.x/iResolution.y;
  vec2 distVec = uv - center;
  distVec.x *= aspectRatio;
  distVec.y *= aspectRatio;
  float distance = length(distVec);
  
  float multiplier = (distance < 1.0) ? ((distance-1.0)*(distance-1.0)) : 0.0;
  float addend = (sin(frequency * distance - t) + centerLight) 
                  * waveStrength * multiplier;
  
  //vec2 newTexCoord = uv + addend*oblique;    
  
  vec4 colorToAdd = sunlightColor * sunlightStrength * addend;
  
  vec3 innercol, excol;
  innercol += (0.05) / rho;
  excol  += (-0.01) / rho;
  return  vec4(innercol, 1.) + colorToAdd ;
}


void main()
{
  
  vec2 uv = ( gl_FragCoord.xy * 2.- iResolution.xy ) / min(iResolution.y,iResolution.x) ;
  float rho = sqrt(uv.x*uv.x + uv.y*uv.y);

  MilkyWay = make_MilkyWay(uv, rho);
  BlackHole = MakeBlackHole(uv, rho);
  
  vec3 from = vec3(1.+iTime*.004,3.+iTime*.002, .5);
  
  vec4 stars = vec4(kali_stars(from, vec3(uv,1.)), 1.); 
 
  float densFog = .95;
  float fogFactor = exp(-(densFog*abs(rho)));
  
  
  gl_FragColor =  stars + fogFactor * (MilkyWay) + BlackHole ;

}
