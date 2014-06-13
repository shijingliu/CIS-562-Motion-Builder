
//////////////////////////////////////////////////////////////////////////
// Transformation.cpp -- Source file for useful Classes about 3D transformations
//
// Liming Zhao
// 11/02/2007
// University of Pennsylvania

#include "Transformation.h"
#define pi 3.1415926
/****************************************************************
*																*
*		    vec2 Member functions								*
*																*
****************************************************************/

// CONSTRUCTORS

vec2::vec2() 
{
}

vec2::vec2(const float x, const float y)
{
        n[VX] = x; n[VY] = y;
}

vec2::vec2(const vec2& v)
{ 
        n[VX] = v.n[VX]; n[VY] = v.n[VY];
}

// ASSIGNMENT OPERATORS

vec2& vec2::operator = (const vec2& v)
{ 
        n[VX] = v.n[VX]; n[VY] = v.n[VY]; return *this;
}

vec2& vec2::operator += ( const vec2& v )
{ 
        n[VX] += v.n[VX]; n[VY] += v.n[VY]; return *this;
}

vec2& vec2::operator -= ( const vec2& v )
{ 
        n[VX] -= v.n[VX]; n[VY] -= v.n[VY]; return *this;
}

vec2& vec2::operator *= ( const float d )
{ 
        n[VX] *= d; n[VY] *= d; return *this;
}

vec2& vec2::operator /= ( const float d )
{ 
        float d_inv = 1.0f/d; n[VX] *= d_inv; n[VY] *= d_inv; return *this;
}

float& vec2::operator [] ( int i) 
{
        assert(!(i < VX || i > VY));		// subscript check
	return n[i];
}

float vec2::operator [] ( int i) const 
{
        assert(!(i < VX || i > VY));
	return n[i];
}


// SPECIAL FUNCTIONS

float vec2::Length() const
{ 
	return sqrt(SqrLength()); 
}

float vec2::SqrLength() const
{ 
        return n[VX]*n[VX] + n[VY]*n[VY];
}

vec2& vec2::Normalize() // it is up to caller to avoid divide-by-zero
{ 
	*this /= Length(); return *this; 
}

// FRIENDS

vec2 operator - (const vec2& a)
{ 
        return vec2(-a.n[VX],-a.n[VY]);
}

vec2 operator + (const vec2& a, const vec2& b)
{ 
        return vec2(a.n[VX]+ b.n[VX], a.n[VY] + b.n[VY]);
}

vec2 operator - (const vec2& a, const vec2& b)
{ 
        return vec2(a.n[VX]-b.n[VX], a.n[VY]-b.n[VY]);
}

vec2 operator * (const vec2& a, const float d)
{ 
        return vec2(d*a.n[VX], d*a.n[VY]);
}

vec2 operator * (const float d, const vec2& a)
{ 
	return a*d; 
}

float operator * (const vec2& a, const vec2& b)
{ 
        return (a.n[VX]*b.n[VX] + a.n[VY]*b.n[VY]);
}

vec2 operator / (const vec2& a, const float d)
{ 
        float d_inv = 1.0f/d; return vec2(a.n[VX]*d_inv, a.n[VY]*d_inv);
}

vec3 operator ^ (const vec2& a, const vec2& b)
{ 
        return vec3(0.0, 0.0, a.n[VX] * b.n[VY] - b.n[VX] * a.n[VY]);
}

int operator == (const vec2& a, const vec2& b)
{ 
        return (a.n[VX] == b.n[VX]) && (a.n[VY] == b.n[VY]);
}

int operator != (const vec2& a, const vec2& b)
{ 
	return !(a == b); 
}

vec2 Prod(const vec2& a, const vec2& b)
{ 
        return vec2(a.n[VX] * b.n[VX], a.n[VY] * b.n[VY]);
}

float Dot(const vec2& a, const vec2& b)
{
	return a*b;
}


/****************************************************************
*																*
*		    vec3 Member functions								*
*																*
****************************************************************/

// CONSTRUCTORS

vec3::vec3() 
{
}

vec3::vec3(const float x, const float y, const float z)
{ 
        n[VX] = x; n[VY] = y; n[VZ] = z;
}

vec3::vec3(const vec3& v)
{ 
        n[VX] = v.n[VX]; n[VY] = v.n[VY]; n[VZ] = v.n[VZ];
}

// ASSIGNMENT OPERATORS

vec3& vec3::operator = (const vec3& v)
{ 
        n[VX] = v.n[VX]; n[VY] = v.n[VY]; n[VZ] = v.n[VZ]; return *this;
}

vec3& vec3::operator += ( const vec3& v )
{ 
        n[VX] += v.n[VX]; n[VY] += v.n[VY]; n[VZ] += v.n[VZ]; return *this;
}

vec3& vec3::operator -= ( const vec3& v )
{ 
        n[VX] -= v.n[VX]; n[VY] -= v.n[VY]; n[VZ] -= v.n[VZ]; return *this;
}

vec3& vec3::operator *= ( const float d )
{ 
        n[VX] *= d; n[VY] *= d; n[VZ] *= d; return *this;
}

vec3& vec3::operator /= ( const float d )
{ 
        float d_inv = 1.0f/d; n[VX] *= d_inv; n[VY] *= d_inv; n[VZ] *= d_inv;
	return *this; 
}

float& vec3::operator [] ( int i) {
        assert(! (i < VX || i > VZ));
	return n[i];
}

float vec3::operator [] ( int i) const {
        assert(! (i < VX || i > VZ));
	return n[i];
}


// SPECIAL FUNCTIONS

float vec3::Length() const
{  
	return sqrt(SqrLength()); 
}

float vec3::SqrLength() const
{  
        return n[VX]*n[VX] + n[VY]*n[VY] + n[VZ]*n[VZ];
}

vec3& vec3::Normalize() // it is up to caller to avoid divide-by-zero
{ 
	*this /= Length(); return *this; 
}

vec3 vec3::Cross(vec3 &v) const
{
	vec3 tmp;
	tmp[0] = n[1] * v.n[2] - n[2] * v.n[1];
	tmp[1] = n[2] * v.n[0] - n[0] * v.n[2];
	tmp[2] = n[0] * v.n[1] - n[1] * v.n[0];
	return tmp;
}

// FRIENDS

vec3 operator - (const vec3& a)
{  
        return vec3(-a.n[VX],-a.n[VY],-a.n[VZ]);
}

vec3 operator + (const vec3& a, const vec3& b)
{ 
        return vec3(a.n[VX]+ b.n[VX], a.n[VY] + b.n[VY], a.n[VZ] + b.n[VZ]);
}

vec3 operator - (const vec3& a, const vec3& b)
{ 
        return vec3(a.n[VX]-b.n[VX], a.n[VY]-b.n[VY], a.n[VZ]-b.n[VZ]);
}

vec3 operator * (const vec3& a, const float d)
{ 
        return vec3(d*a.n[VX], d*a.n[VY], d*a.n[VZ]);
}

vec3 operator * (const float d, const vec3& a)
{ 
	return a*d; 
}

vec3 operator * (const mat3& a, const vec3& v) 
{
#define ROWCOL(i) a.v[i].n[0]*v.n[VX] + a.v[i].n[1]*v.n[VY] \
        + a.v[i].n[2]*v.n[VZ]
	return vec3(ROWCOL(0), ROWCOL(1), ROWCOL(2));
#undef ROWCOL
}

float operator * (const vec3& a, const vec3& b)
{ 
        return (a.n[VX]*b.n[VX] + a.n[VY]*b.n[VY] + a.n[VZ]*b.n[VZ]);
}

vec3 operator / (const vec3& a, const float d)
{ 
	float d_inv = 1.0f/d; 
        return vec3(a.n[VX]*d_inv, a.n[VY]*d_inv, a.n[VZ]*d_inv);
}

vec3 operator ^ (const vec3& a, const vec3& b) 
{
        return vec3(a.n[VY]*b.n[VZ] - a.n[VZ]*b.n[VY],
                a.n[VZ]*b.n[VX] - a.n[VX]*b.n[VZ],
                a.n[VX]*b.n[VY] - a.n[VY]*b.n[VX]);
}

int operator == (const vec3& a, const vec3& b)
{ 
        return (a.n[VX] == b.n[VX]) && (a.n[VY] == b.n[VY]) && (a.n[VZ] == b.n[VZ]);
}

int operator != (const vec3& a, const vec3& b)
{ 
	return !(a == b); 
}

vec3 Prod(const vec3& a, const vec3& b)
{ 
        return vec3(a.n[VX] * b.n[VX], a.n[VY] * b.n[VY], a.n[VZ] * b.n[VZ]);
}

float Dot(const vec3& a, const vec3& b)
{
	return a*b;
}
/****************************************************************
*																*
*		    vec4 Member functions								*
*																*
****************************************************************/

// CONSTRUCTORS

vec4::vec4()
{
}

vec4::vec4(const Real x, const Real y, const Real z, const Real w)
{
        n[VX] = x; n[VY] = y; n[VZ] = z; n[VW] = w;
}

vec4::vec4(const Real d)
{
        n[VX] = n[VY] = n[VZ] = n[VW] = d;
}

vec4::vec4(const vec4& v)
{
        n[VX] = v.n[VX]; n[VY] = v.n[VY]; n[VZ] = v.n[VZ]; n[VW] = v.n[VW];
}

vec4::vec4(const vec3& v)
{
        n[VX] = v[VX]; n[VY] = v[VY]; n[VZ] = v[VZ]; n[VW] = 1.0;
}

vec4::vec4(const vec3& v, const Real d)
{
        n[VX] = v[VX]; n[VY] = v[VY]; n[VZ] = v[VZ];  n[VW] = d;
}


// ASSIGNMENT OPERATORS

vec4& vec4::operator = (const vec4& v)
{
        n[VX] = v.n[VX]; n[VY] = v.n[VY]; n[VZ] = v.n[VZ]; n[VW] = v.n[VW];
        return *this;
}

vec4& vec4::operator += ( const vec4& v )
{
        n[VX] += v.n[VX]; n[VY] += v.n[VY]; n[VZ] += v.n[VZ]; n[VW] += v.n[VW];
        return *this;
}

vec4& vec4::operator -= ( const vec4& v )
{
        n[VX] -= v.n[VX]; n[VY] -= v.n[VY]; n[VZ] -= v.n[VZ]; n[VW] -= v.n[VW];
        return *this;
}

vec4& vec4::operator *= ( const Real d )
{
        n[VX] *= d; n[VY] *= d; n[VZ] *= d; n[VW] *= d;
        return *this;
}

vec4& vec4::operator /= ( const Real d )
{
        Real d_inv = 1./d;
        n[VX] *= d_inv; n[VY] *= d_inv; n[VZ] *= d_inv; n[VW] *= d_inv;
        return *this;
}

Real& vec4::operator [] ( int i)
{
        assert(! (i < VX || i > VW));
        return n[i];
}

Real vec4::operator [] ( int i) const
{
        assert(! (i < VX || i > VW));
        return n[i];
}

// SPECIAL FUNCTIONS

Real vec4::Length() const
{
        return sqrt(Length2());
}

Real vec4::Length2() const
{
        return n[VX]*n[VX] + n[VY]*n[VY] + n[VZ]*n[VZ] + n[VW]*n[VW];
}

vec4& vec4::Normalize() // it is up to caller to avoid divide-by-zero
{
        *this /= Length(); return *this;
}

int vec4::dim() const								// SHL added - returns dimension of vector
{
        return (sizeof(n)/sizeof(Real));
}

// FRIENDS

vec4 operator - (const vec4& a)
{
        return vec4(-a.n[VX],-a.n[VY],-a.n[VZ],-a.n[VW]);
}

vec4 operator + (const vec4& a, const vec4& b)
{
        return vec4(a.n[VX] + b.n[VX], a.n[VY] + b.n[VY], a.n[VZ] + b.n[VZ], a.n[VW] + b.n[VW]);
}

vec4 operator - (const vec4& a, const vec4& b)
{
        return vec4(a.n[VX] - b.n[VX], a.n[VY] - b.n[VY], a.n[VZ] - b.n[VZ], a.n[VW] - b.n[VW]);
}

vec4 operator * (const vec4& a, const Real d)
{
        return vec4(d*a.n[VX], d*a.n[VY], d*a.n[VZ], d*a.n[VW] );
}

vec4 operator * (const Real d, const vec4& a)
{
        return a*d;
}

vec4 operator * (const mat4& a, const vec4& v)
{
#define ROWCOL(i) a.v[i].n[0]*v.n[VX] + a.v[i].n[1]*v.n[VY] \
        + a.v[i].n[2]*v.n[VZ] + a.v[i].n[3]*v.n[VW]
        return vec4(ROWCOL(0), ROWCOL(1), ROWCOL(2), ROWCOL(3));
#undef ROWCOL // (i)
}

vec4 operator * (const vec4& v, const mat4& a)
{
        return a.transpose() * v;
}

Real operator * (const vec4& a, const vec4& b)
{
        return (a.n[VX]*b.n[VX] + a.n[VY]*b.n[VY] + a.n[VZ]*b.n[VZ] + a.n[VW]*b.n[VW]);
}

vec4 operator / (const vec4& a, const Real d)
{
        Real d_inv = 1./d;
        return vec4(a.n[VX]*d_inv, a.n[VY]*d_inv, a.n[VZ]*d_inv, a.n[VW]*d_inv);
}

int operator == (const vec4& a, const vec4& b)
{
        return (a.n[VX] == b.n[VX]) && (a.n[VY] == b.n[VY]) && (a.n[VZ] == b.n[VZ]) && (a.n[VW] == b.n[VW]);
}

int operator != (const vec4& a, const vec4& b)
{
        return !(a == b);
}

#ifdef ALGEBRAIOSTREAMS
ostream& operator << (ostream& s, const vec4& v)
{
        return s << "[ " << v.n[VX] << ' ' << v.n[VY] << ' ' << v.n[VZ] << ' ' << v.n[VW] << " ]";
}

// stream& operator >> (istream& s, vec4& v)
//{
//    vec4	v_tmp;
//    char	c = ' ';
//
//    while (isspace(c))
//	s >> c;
//    // The vectors can be formatted either as x y z w or | x y z w |
//    if (c == '[') {
//	s >> v_tmp[VX] >> v_tmp[VY] >> v_tmp[VZ] >> v_tmp[VW];
//	while (s >> c && isspace(c)) ;
//	if (c != ']')
//	    s.set(_bad);
//	}
//    else {
//	s.putback(c);
//	s >> v_tmp[VX] >> v_tmp[VY] >> v_tmp[VZ] >> v_tmp[VW];
//	}
//    if (s)
//	v = v_tmp;
//    return s;
//}
#endif // ALGEBRAIOSTREAMS

void swap(vec4& a, vec4& b)
{
        vec4 tmp(a); a = b; b = tmp;
}

vec4 min(const vec4& a, const vec4& b)
{
        return vec4(MIN(a.n[VX], b.n[VX]), MIN(a.n[VY], b.n[VY]), MIN(a.n[VZ], b.n[VZ]), MIN(a.n[VW], b.n[VW]));
}

vec4 max(const vec4& a, const vec4& b)
{
        return vec4(MAX(a.n[VX], b.n[VX]), MAX(a.n[VY], b.n[VY]), MAX(a.n[VZ], b.n[VZ]), MAX(a.n[VW], b.n[VW]));
}

vec4 prod(const vec4& a, const vec4& b)
{
        return vec4(a.n[VX] * b.n[VX], a.n[VY] * b.n[VY], a.n[VZ] * b.n[VZ], a.n[VW] * b.n[VW]);
}

/****************************************************************
*																*
*		    mat3 member functions								*
*																*
****************************************************************/  

// CONSTRUCTORS

mat3::mat3() 
{
	v[0] = vec3(0.0,0.0,0.0);
	v[1] = v[2] = v[0];
}

mat3::mat3(const vec3& v0, const vec3& v1, const vec3& v2)
{ 
	v[0] = v0; v[1] = v1; v[2] = v2; 
}

mat3::mat3(const mat3& m)
{ 
	v[0] = m.v[0]; v[1] = m.v[1]; v[2] = m.v[2]; 
}

// Static functions

mat3 mat3::Identity()
{
	return mat3(vec3(1.0, 0.0, 0.0),
		vec3(0.0, 1.0, 0.0),
		vec3(0.0, 0.0, 1.0));
}

mat3 mat3::Rotation3DRad(const vec3& axis, const float angleRad)
{
	float c = cos(angleRad), s = sin(angleRad), t = 1.0f - c;
	vec3 Axis = axis;
	Axis.Normalize();
        return mat3(vec3(t * Axis[VX] * Axis[VX] + c,
                t * Axis[VX] * Axis[VY] - s * Axis[VZ],
                t * Axis[VX] * Axis[VZ] + s * Axis[VY]),
                vec3(t * Axis[VX] * Axis[VY] + s * Axis[VZ],
                t * Axis[VY] * Axis[VY] + c,
                t * Axis[VY] * Axis[VZ] - s * Axis[VX]),
                vec3(t * Axis[VX] * Axis[VZ] - s * Axis[VY],
                t * Axis[VY] * Axis[VZ] + s * Axis[VX],
                t * Axis[VZ] * Axis[VZ] + c)
		);
}

// ASSIGNMENT OPERATORS

mat3& mat3::operator = ( const mat3& m )
{ 
	v[0] = m.v[0]; v[1] = m.v[1]; v[2] = m.v[2]; 
	return *this; 
}

mat3& mat3::operator += ( const mat3& m )
{ 
	v[0] += m.v[0]; v[1] += m.v[1]; v[2] += m.v[2]; 
	return *this; 
}

mat3& mat3::operator -= ( const mat3& m )
{ 
	v[0] -= m.v[0]; v[1] -= m.v[1]; v[2] -= m.v[2]; 
	return *this; 
}

mat3& mat3::operator *= ( const float d )
{ 
	v[0] *= d; v[1] *= d; v[2] *= d; 
	return *this; 
}

mat3& mat3::operator /= ( const float d )
{ 
	v[0] /= d; v[1] /= d; v[2] /= d; 
	return *this; 
}

vec3& mat3::operator [] ( int i) 
{
        assert(! (i < VX || i > VZ));
	return v[i];
}

const vec3& mat3::operator [] ( int i) const 
{
        assert(!(i < VX || i > VZ));
	return v[i];
}

// SPECIAL FUNCTIONS

mat3 mat3::Transpose() const 
{
	return mat3(vec3(v[0][0], v[1][0], v[2][0]),
		vec3(v[0][1], v[1][1], v[2][1]),
		vec3(v[0][2], v[1][2], v[2][2]));
}

// FRIENDS

mat3 operator - (const mat3& a)
{ 
	return mat3(-a.v[0], -a.v[1], -a.v[2]); 
}

mat3 operator + (const mat3& a, const mat3& b)
{ 
	return mat3(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]); 
}

mat3 operator - (const mat3& a, const mat3& b)
{ 
	return mat3(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]); 
}

mat3 operator * (const mat3& a, const mat3& b) 
{
#define ROWCOL(i, j) \
	a.v[i].n[0]*b.v[0][j] + a.v[i].n[1]*b.v[1][j] + a.v[i].n[2]*b.v[2][j]
	return mat3(vec3(ROWCOL(0,0), ROWCOL(0,1), ROWCOL(0,2)),
		vec3(ROWCOL(1,0), ROWCOL(1,1), ROWCOL(1,2)),
		vec3(ROWCOL(2,0), ROWCOL(2,1), ROWCOL(2,2)));
#undef ROWCOL // (i, j)
}

mat3 operator * (const mat3& a, const float d)
{ 
	return mat3(a.v[0] * d, a.v[1] * d, a.v[2] * d); 
}

mat3 operator * (const float d, const mat3& a)
{ 
	return a*d; 
}

mat3 operator / (const mat3& a, const float d)
{ 
	return mat3(a.v[0] / d, a.v[1] / d, a.v[2] / d); 
}

int operator == (const mat3& a, const mat3& b)
{ 
	return (a.v[0] == b.v[0]) && (a.v[1] == b.v[1]) && (a.v[2] == b.v[2]); 
}

int operator != (const mat3& a, const mat3& b)
{ 
	return !(a == b); 
}
/****************************************************************
*																*
*		    mat4 member functions								*
*																*
****************************************************************/

// CONSTRUCTORS

mat4::mat4()
{
}

mat4::mat4(const vec4& v0, const vec4& v1, const vec4& v2, const vec4& v3)
{
        v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
}

mat4::mat4(const Real d)
{
        v[0] = v[1] = v[2] = v[3] = vec4(d);
}

mat4::mat4(const mat4& m)
{
        v[0] = m.v[0]; v[1] = m.v[1]; v[2] = m.v[2]; v[3] = m.v[3];
}

mat4::mat4(const Real* d)
{
        v[0] = vec4(d[0], d[4], d[8], d[12]);
        v[1] = vec4(d[1], d[5], d[9], d[13]);
        v[2] = vec4(d[2], d[6], d[10], d[14]);
        v[3] = vec4(d[3], d[7], d[11], d[15]);
}

mat4::mat4(const mat3& m)
{
        v[0] = vec4(m[0], 0);
        v[1] = vec4(m[1], 0);
        v[2] = vec4(m[2], 0);
        v[3] = vec4(0, 0, 0, 1);
}

mat4::mat4(const mat3& m, const vec3& t)
{
        v[0] = vec4(m[0], t[0]);
        v[1] = vec4(m[1], t[1]);
        v[2] = vec4(m[2], t[2]);
        v[3] = vec4(0, 0, 0, 1);
}

// Static functions

mat4 mat4::identity()
{
        return mat4(vec4(1.0, 0.0, 0.0, 0.0),
                vec4(0.0, 1.0, 0.0, 0.0),
                vec4(0.0, 0.0, 1.0, 0.0),
                vec4(0.0, 0.0, 0.0, 1.0));
}

mat4 mat4::translation3D(const vec3& v)
{
        return mat4(vec4(1.0, 0.0, 0.0, v[VX]),
                vec4(0.0, 1.0, 0.0, v[VY]),
                vec4(0.0, 0.0, 1.0, v[VZ]),
                vec4(0.0, 0.0, 0.0, 1.0));
}

mat4 mat4::rotation3DDeg(const vec3& axis, const Real angleDeg)
{
        Real angleRad = angleDeg * Deg2Rad;
        return rotation3DRad(axis, angleRad);
}

mat4 mat4::rotation3DRad(const vec3& axis, const Real angleRad)
{
        Real  c = cos(angleRad),
                s = sin(angleRad),
                t = 1.0 - c;
        vec3 Axis = axis;
        Axis.Normalize();
        return mat4(vec4(t * Axis[VX] * Axis[VX] + c,
                t * Axis[VX] * Axis[VY] - s * Axis[VZ],
                t * Axis[VX] * Axis[VZ] + s * Axis[VY],
                0.0),
                vec4(t * Axis[VX] * Axis[VY] + s * Axis[VZ],
                t * Axis[VY] * Axis[VY] + c,
                t * Axis[VY] * Axis[VZ] - s * Axis[VX],
                0.0),
                vec4(t * Axis[VX] * Axis[VZ] - s * Axis[VY],
                t * Axis[VY] * Axis[VZ] + s * Axis[VX],
                t * Axis[VZ] * Axis[VZ] + c,
                0.0),
                vec4(0.0, 0.0, 0.0, 1.0));
}

mat4 mat4::scaling3D(const vec3& scaleVector)
{
        return mat4(vec4(scaleVector[VX], 0.0, 0.0, 0.0),
                vec4(0.0, scaleVector[VY], 0.0, 0.0),
                vec4(0.0, 0.0, scaleVector[VZ], 0.0),
                vec4(0.0, 0.0, 0.0, 1.0));
}

mat4 mat4::perspective3D(const Real d)
{
        return mat4(vec4(1.0, 0.0, 0.0, 0.0),
                vec4(0.0, 1.0, 0.0, 0.0),
                vec4(0.0, 0.0, 1.0, 0.0),
                vec4(0.0, 0.0, 1.0/d, 0.0));
}

// ASSIGNMENT OPERATORS

mat4& mat4::operator = ( const mat4& m )
{
        v[0] = m.v[0]; v[1] = m.v[1]; v[2] = m.v[2]; v[3] = m.v[3];
        return *this;
}

mat4& mat4::operator += ( const mat4& m )
{
        v[0] += m.v[0]; v[1] += m.v[1]; v[2] += m.v[2]; v[3] += m.v[3];
        return *this;
}

mat4& mat4::operator -= ( const mat4& m )
{
        v[0] -= m.v[0]; v[1] -= m.v[1]; v[2] -= m.v[2]; v[3] -= m.v[3];
        return *this;
}

mat4& mat4::operator *= ( const Real d )
{
        v[0] *= d; v[1] *= d; v[2] *= d; v[3] *= d;
        return *this;
}

mat4& mat4::operator /= ( const Real d )
{
        v[0] /= d; v[1] /= d; v[2] /= d; v[3] /= d;
        return *this;
}

vec4& mat4::operator [] ( int i)
{
        assert(! (i < VX || i > VW));
        return v[i];
}

const vec4& mat4::operator [] ( int i) const
{
        assert(! (i < VX || i > VW));
        return v[i];
}

// SPECIAL FUNCTIONS;

mat4 mat4::transpose() const
{
        return mat4(vec4(v[0][0], v[1][0], v[2][0], v[3][0]),
                vec4(v[0][1], v[1][1], v[2][1], v[3][1]),
                vec4(v[0][2], v[1][2], v[2][2], v[3][2]),
                vec4(v[0][3], v[1][3], v[2][3], v[3][3]));
}

void mat4::getData(Real* d)
{
        d[0] = v[0][0]; d[1] = v[1][0]; d[2] = v[2][0]; d[3] = v[3][0];
        d[4] = v[0][1]; d[5] = v[1][1]; d[6] = v[2][1]; d[7] = v[3][1];
        d[8] = v[0][2]; d[9] = v[1][2]; d[10] = v[2][2]; d[11] = v[3][2];
        d[12] = v[0][3]; d[13] = v[1][3]; d[14] = v[2][3]; d[15] = v[3][3];
}

// FRIENDS

mat4 operator - (const mat4& a)
{
        return mat4(-a.v[0], -a.v[1], -a.v[2], -a.v[3]);
}

mat4 operator + (const mat4& a, const mat4& b)
{
        return mat4(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2], a.v[3] + b.v[3]);
}

mat4 operator - (const mat4& a, const mat4& b)
{
        return mat4(a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2], a.v[3] - b.v[3]);
}

mat4 operator * (const mat4& a, const mat4& b)
{
#define ROWCOL(i, j) a.v[i].n[0]*b.v[0][j] + a.v[i].n[1]*b.v[1][j] + \
        a.v[i].n[2]*b.v[2][j] + a.v[i].n[3]*b.v[3][j]
        return mat4(
                vec4(ROWCOL(0,0), ROWCOL(0,1), ROWCOL(0,2), ROWCOL(0,3)),
                vec4(ROWCOL(1,0), ROWCOL(1,1), ROWCOL(1,2), ROWCOL(1,3)),
                vec4(ROWCOL(2,0), ROWCOL(2,1), ROWCOL(2,2), ROWCOL(2,3)),
                vec4(ROWCOL(3,0), ROWCOL(3,1), ROWCOL(3,2), ROWCOL(3,3))
                );
#undef ROWCOL
}

mat4 operator * (const mat4& a, const Real d)
{
        return mat4(a.v[0] * d, a.v[1] * d, a.v[2] * d, a.v[3] * d);
}

mat4 operator * (const Real d, const mat4& a)
{
        return a*d;
}

mat4 operator / (const mat4& a, const Real d)
{
        return mat4(a.v[0] / d, a.v[1] / d, a.v[2] / d, a.v[3] / d);
}

int operator == (const mat4& a, const mat4& b)
{
        return ((a.v[0] == b.v[0]) && (a.v[1] == b.v[1]) && (a.v[2] == b.v[2]) && (a.v[3] == b.v[3]));
}

int operator != (const mat4& a, const mat4& b)
{
        return !(a == b);
}

#ifdef ALGEBRAIOSTREAMS
ostream& operator << (ostream& s, const mat4& m)
{
        return s << m.v[VX] << '\n' << m.v[VY] << '\n' << m.v[VZ] << '\n' << m.v[VW];
}

// istream& operator >> (istream& s, mat4& m)
//{
//	mat4 m_tmp;
//	s >> m_tmp[VX] >> m_tmp[VY] >> m_tmp[VZ] >> m_tmp[VW];
//	if (s)
//		m = m_tmp;
//	return s;
//}
#endif // ALGEBRAIOSTREAMS

void swap(mat4& a, mat4& b)
{
        mat4 tmp(a); a = b; b = tmp;
}

mat4 mat4::inverse()	const    // Gauss-Jordan elimination with partial pivoting
{
        mat4 a(*this),	    // As a evolves from original mat into identity
                b(mat4::identity());   // b evolves from identity into inverse(a)
        int i, j, i1;

        // Loop over cols of a from left to right, eliminating above and below diag
        for (j=0; j<4; j++) {   // Find largest pivot in column j among rows j..3
                i1 = j;		    // Row with largest pivot candidate
                for (i=j+1; i<4; i++)
                        if (fabs(a.v[i].n[j]) > fabs(a.v[i1].n[j]))
                                i1 = i;

                // Swap rows i1 and j in a and b to put pivot on diagonal
                swap(a.v[i1], a.v[j]);
                swap(b.v[i1], b.v[j]);

                // Scale row j to have a unit diagonal
                if (a.v[j].n[j]==0.)
                        ALGEBRA_ERROR("mat4::inverse: singular matrix; can't invert\n");
                b.v[j] /= a.v[j].n[j];
                a.v[j] /= a.v[j].n[j];

                // Eliminate off-diagonal elems in col j of a, doing identical ops to b
                for (i=0; i<4; i++)
                        if (i!=j) {
                                b.v[i] -= a.v[i].n[j]*b.v[j];
                                a.v[i] -= a.v[i].n[j]*a.v[j];
                        }
        }
        return b;
}

/****************************************************************
*																*
*		    Quaternion member functions							*
*																*
****************************************************************/

// CONSTRUCTORS

Quaternion::Quaternion()
{
}

Quaternion::Quaternion(const float w, const float x, const float y, const float z)
{
        n[VW] = w; n[VX] = x; n[VY] = y; n[VZ] = z;
}

Quaternion::Quaternion(const Quaternion& q)
{
        n[VW] = q.n[VW]; n[VX] = q.n[VX]; n[VY] = q.n[VY]; n[VZ] = q.n[VZ];
}

// Static functions

float Quaternion::Dot(const Quaternion& q0, const Quaternion& q1)
{
        return q0.n[VW] * q1.n[VW] + q0.n[VX] * q1.n[VX] + q0.n[VY] * q1.n[VY] + q0.n[VZ] * q1.n[VZ];
}

Quaternion Quaternion::UnitInverse(const Quaternion& q)
{
        return Quaternion(q.n[VW], -q.n[VX], -q.n[VY], -q.n[VZ]);
}

float Quaternion::CounterWarp(float t, float fCos)
{
	const float ATTENUATION = 0.82279687f;
	const float WORST_CASE_SLOPE = 0.58549219f;

	float fFactor = 1.0f - ATTENUATION * fCos;
	fFactor *= fFactor;
	float fK = WORST_CASE_SLOPE * fFactor;

	return t * (fK * t * (2.0f * t - 3.0f) + 1.0f + fK);
}

static const float ISQRT_NEIGHBORHOOD = 0.959066f;
static const float ISQRT_SCALE = 1.000311f;
static const float ISQRT_ADDITIVE_CONSTANT = ISQRT_SCALE / (float)sqrt(ISQRT_NEIGHBORHOOD);
static const float ISQRT_FACTOR = ISQRT_SCALE * (-0.5f / (ISQRT_NEIGHBORHOOD * (float)sqrt(ISQRT_NEIGHBORHOOD)));
float Quaternion::ISqrt_approx_in_neighborhood(float s)
{
	return ISQRT_ADDITIVE_CONSTANT + (s - ISQRT_NEIGHBORHOOD) * ISQRT_FACTOR;	
}

// Assignment operators

Quaternion& Quaternion::operator = (const Quaternion& q)
{
        n[VW] = q.n[VW]; n[VX] = q.n[VX]; n[VY] = q.n[VY]; n[VZ] = q.n[VZ];
	return *this;
}

Quaternion& Quaternion::operator += (const Quaternion& q)
{
        n[VW] += q.n[VW]; n[VX] += q.n[VX]; n[VY] += q.n[VY]; n[VZ] += q.n[VZ];
	return *this;
}

Quaternion& Quaternion::operator -= (const Quaternion& q)
{
        n[VW] -= q.n[VW]; n[VX] -= q.n[VX]; n[VY] -= q.n[VY]; n[VZ] -= q.n[VZ];
	return *this;
}

Quaternion& Quaternion::operator *= (const Quaternion& q)
{
        *this = Quaternion(n[VW] * q.n[VW] - n[VX] * q.n[VX] - n[VY] * q.n[VY] - n[VZ] * q.n[VZ],
                n[VW] * q.n[VX] + n[VX] * q.n[VW] + n[VY] * q.n[VZ] - n[VZ] * q.n[VY],
                n[VW] * q.n[VY] + n[VY] * q.n[VW] + n[VZ] * q.n[VX] - n[VX] * q.n[VZ],
                n[VW] * q.n[VZ] + n[VZ] * q.n[VW] + n[VX] * q.n[VY] - n[VY] * q.n[VX]);
	return *this;
}

Quaternion& Quaternion::operator *= (const float d)
{
        n[VW] *= d; n[VX] *= d;	n[VY] *= d; n[VZ] *= d;
	return *this;
}

Quaternion& Quaternion::operator /= (const float d)
{
        n[VW] /= d; n[VX] /= d;	n[VY] /= d; n[VZ] /= d;
	return *this;
}

// Indexing
float& Quaternion::operator [](int i)
{
	return n[i];
}

float Quaternion::operator [](int i) const
{
	return n[i];
}

float& Quaternion::W()
{
        return n[VW];
}

float Quaternion::W() const
{
        return n[VW];
}

float& Quaternion::X()
{
        return n[VX];
}

float Quaternion::X() const
{
        return n[VX];
}

float& Quaternion::Y()
{
        return n[VY];
}

float Quaternion::Y() const
{
        return n[VY];
}

float& Quaternion::Z()
{
        return n[VZ];
}

float Quaternion::Z() const
{
        return n[VZ];
}

// Friends

Quaternion operator - (const Quaternion& q)
{
        return Quaternion(-q.n[VW], -q.n[VX], -q.n[VY], -q.n[VZ]);
}

Quaternion operator + (const Quaternion& q0, Quaternion& q1)
{
        return Quaternion(q0.n[VW] + q1.n[VW], q0.n[VX] + q1.n[VX], q0.n[VY] + q1.n[VY], q0.n[VZ] + q1.n[VZ]);
}

Quaternion operator - (const Quaternion& q0, const Quaternion& q1)
{
        return Quaternion(q0.n[VW] - q1.n[VW], q0.n[VX] - q1.n[VX], q0.n[VY] - q1.n[VY], q0.n[VZ] - q1.n[VZ]);
}

Quaternion operator * (const Quaternion& q, const float d)
{
        return Quaternion(q.n[VW] * d, q.n[VX] * d, q.n[VY] * d, q.n[VZ] * d);
}

Quaternion operator * (const float d, const Quaternion& q)
{
        return Quaternion(q.n[VW] * d, q.n[VX] * d, q.n[VY] * d, q.n[VZ] * d);
}

Quaternion operator * (const Quaternion& q0, const Quaternion& q1)
{
        return Quaternion(q0.n[VW] * q1.n[VW] - q0.n[VX] * q1.n[VX] - q0.n[VY] * q1.n[VY] - q0.n[VZ] * q1.n[VZ],
                q0.n[VW] * q1.n[VX] + q0.n[VX] * q1.n[VW] + q0.n[VY] * q1.n[VZ] - q0.n[VZ] * q1.n[VY],
                q0.n[VW] * q1.n[VY] + q0.n[VY] * q1.n[VW] + q0.n[VZ] * q1.n[VX] - q0.n[VX] * q1.n[VZ],
                q0.n[VW] * q1.n[VZ] + q0.n[VZ] * q1.n[VW] + q0.n[VX] * q1.n[VY] - q0.n[VY] * q1.n[VX]);
}

Quaternion operator / (const Quaternion& q, const float d)
{
        return Quaternion(q.n[VW] / d, q.n[VX] / d, q.n[VY] / d, q.n[VZ] / d);
}

bool operator == (const Quaternion& q0, const Quaternion& q1)
{
        return (q0.n[VW] == q1.n[VW]) && (q0.n[VX] == q1.n[VX]) && (q0.n[VY] == q1.n[VY]) && (q0.n[VZ] == q1.n[VZ]);
}

bool operator != (const Quaternion& q0, const Quaternion& q1)
{
	return !(q0 == q1); 
}

// special functions

float Quaternion::SqrLength() const
{
        return n[VW] * n[VW] + n[VX] * n[VX] + n[VY] * n[VY] + n[VZ] * n[VZ];
}

float Quaternion::Length() const
{
	return sqrt(SqrLength());
}

Quaternion& Quaternion::Normalize()
{
	float l = Length();
	if (l < EPSILON || abs(l) > 1e6)
	{
		FromAxisAngle(axisY, 0.0f);
	}else
	{
		*this /= l;
	}

	return *this; 
}

Quaternion& Quaternion::FastNormalize() 
{
        float s = n[VW] * n[VW] + n[VX] * n[VX] + n[VY] * n[VY] + n[VZ] * n[VZ]; // length^2
	float k = ISqrt_approx_in_neighborhood(s);

	if (s <= 0.91521198) {
		k *= ISqrt_approx_in_neighborhood(k * k * s);

		if (s <= 0.65211970) {
			k *= ISqrt_approx_in_neighborhood(k * k * s);
		}
	}

        n[VW] *= k;
        n[VX] *= k;
        n[VY] *= k;
        n[VZ] *= k;

	return * this;
}

Quaternion Quaternion::Inverse() const
{
        return Quaternion(n[VW], -n[VX], -n[VY], -n[VZ]);
}

Quaternion Quaternion::Exp(const Quaternion& q)
{
	// q = A*(x*i+y*j+z*k) where (x,y,z) is unit length
	// exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k)
        float angle = sqrt(q.n[VX] * q.n[VX] + q.n[VY] * q.n[VY] + q.n[VZ] * q.n[VZ]);
	float sn, cs;
	sn = sin(angle);
	cs = cos(angle);

	// When A is near zero, sin(A)/A is approximately 1.  Use
	// exp(q) = cos(A)+A*(x*i+y*j+z*k)
	float coeff = ( abs(sn) < EPSILON ? 1.0f : sn/angle );

        Quaternion result(cs, coeff * q.n[VX], coeff * q.n[VY], coeff * q.n[VZ]);

	return result;
}

Quaternion Quaternion::Log(const Quaternion& q)
{
	// q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length
	// log(q) = A*(x*i+y*j+z*k)

        float angle = acos(q.n[VW]);
	float sn = sin(angle);

	// When A is near zero, A/sin(A) is approximately 1.  Use
	// log(q) = sin(A)*(x*i+y*j+z*k)
	float coeff = ( abs(sn) < EPSILON ? 1.0f : angle/sn );

        return Quaternion(0.0f, coeff * q.n[VX], coeff * q.n[VY], coeff * q.n[VZ]);
}

void Quaternion::Zero()
{
        n[VW] = n[VX] = n[VY] = n[VZ] = 0.0f;
}

// Conversion functions
void Quaternion::ToAxisAngle (vec3& axis, float& angleRad) const
{
	float fLength = Length();

	if ( fLength < EPSILON )
	{
		angleRad = 0;
                axis[VX] = 0;
                axis[VY] = 0;
                axis[VZ] = 0;
	}
	else
	{
                angleRad = 2.0f * acos(n[VW]);
		float invLength = 1.0f / fLength;
                axis[VX] = n[VX] * invLength;
                axis[VY] = n[VY] * invLength;
                axis[VZ] = n[VZ] * invLength;
	}
}

void Quaternion::FromAxisAngle (const vec3& axis, float angleRad)
{
	float fHalfAngle = angleRad * 0.5f;
	float sn = sin(fHalfAngle);
        n[VW] = cos(fHalfAngle);
        n[VX] = axis[VX] * sn;
        n[VY] = axis[VY] * sn;
        n[VZ] = axis[VZ] * sn;
}

/****************************************************************
*																*
*		    Transform member functions							*
*																*
****************************************************************/   
 
// Constructors
Transform::Transform()
{
	Identity();
}

Transform::Transform(const vec3& translation, const mat3& rotation)
{
	m_translation = translation;
	m_rotation = rotation;
}  

Transform::Transform(const vec3& translation)
{
	m_translation = translation;
	m_rotation = mat3::Identity();
}

Transform::Transform(const mat3& rotation)
{
	m_translation = vec3(0.0f, 0.0f, 0.0f);
	m_rotation = rotation;
}

Transform::Transform(const Transform& transform)
{
	m_translation = transform.m_translation;
	m_rotation = transform.m_rotation;
}

// Destructor
Transform::~Transform(void)
{
}

// Member functions
void Transform::Identity()
{
	m_translation = vec3(0.0f, 0.0f, 0.0f);
	m_rotation = mat3::Identity();
}

Transform& Transform::operator = (const Transform& source)
{
	m_translation = source.m_translation;
	m_rotation = source.m_rotation;
	return *this;
}

Transform Transform::Lerp(const float fPerc, const Transform& t0, const Transform& t1)
{
	Transform result;
	result.m_translation = t0.m_translation * (1.0f - fPerc) + t1.m_translation * fPerc;

	Quaternion q0, q1, q;
	q0.FromRotation(t0.m_rotation);
	q1.FromRotation(t1.m_rotation);
	float d = Quaternion::Dot(q0, q1);
	if (d < 0)
		q1 = -q1;
	q = Quaternion::Slerp(fPerc, q0, q1);
	result.m_rotation = q.ToRotation();

	return result;
}

//************************************************************************
// Functions you need to implement    
//************************************************************************

//======================================================================== 
// 3D Matrix

//////////////////////////////////////////////////////////////////////////
// Convert a rotation matrix to Euler angles in ZXY rotation order
// Euler angles are stored to angleRad in radians
// Return true if there is no Gimbol lock, otherwise false
bool mat3::ToEulerAnglesZXY(vec3& angleRad) const
{
    float theta_x=asin(float(v[2].n[1]));
	float theta_y,theta_z;
	if(cos(theta_x)!=0)
	{
		theta_y=atan2(-v[2].n[0]/cos(theta_x),v[2].n[2]/cos(theta_x));
		theta_z=atan2(-v[0].n[1]/cos(theta_x),v[1].n[1]/cos(theta_x));
	}
	else
	{
		if(theta_x==pi/2.0)
		{
			theta_y=0;
			theta_z=atan2(v[0].n[2],v[0].n[0]);
		}
		else
		{
			theta_z=0;
			theta_y=atan2(v[0].n[2],v[0].n[0]);
		}
	}
	// Replace the following code with your code
	angleRad = vec3(theta_x,theta_y,theta_z);
	return true;
}

//////////////////////////////////////////////////////////////////////////   
// Convert Euler angles in ZXY rotation order to rotation matrix
// The input Euler angles are in angleRad in radians  
// The returning mat3 is the corresponding rotation matrix
mat3 mat3::FromEulerAnglesZXY(const vec3& angleRad)
{
	// Replace the following code with your code
	
	double thetaX = angleRad[0]; 
	double thetaY = angleRad[1]; 
	double thetaZ = angleRad[2];   

	v[0][0] = cos(thetaZ)*cos(thetaY) - sin(thetaZ)*sin(thetaX)*sin(thetaY);  
	v[0][1] = -sin(thetaZ)*cos(thetaX);
	v[0][2] = cos(thetaZ)*sin(thetaY) + sin(thetaZ)*sin(thetaX)*cos(thetaY); 
	v[1][0] = sin(thetaZ)*cos(thetaY) + cos(thetaZ)*sin(thetaX)*sin(thetaY); 
	v[1][1] = cos(thetaZ)*cos(thetaX); 
	v[1][2] = sin(thetaZ)*sin(thetaY) - cos(thetaZ)*sin(thetaX)*cos(thetaY); 
	v[2][0] = -cos(thetaX)*sin(thetaY); 
	v[2][1] = sin(thetaX); 
	v[2][2] = cos(thetaX)*cos(thetaY);

	return *this; 
}

//////////////////////////////////////////////////////////////////////////
// Convert rotation matrix to Quaternion
// Return the corresponding Quaternion
Quaternion mat3::ToQuaternion() const
{
	// Replace the following code with your code
	double S = sqrt(1+v[0][0]+v[1][1]+v[2][2])/2.0; 
    double Vx =  (v[2][1] - v[1][2])/(4*S); 
	double Vy =  (v[0][2] - v[2][0])/(4*S); 
	double Vz =  (v[1][0] - v[0][1])/(4*S); 

	Quaternion q;
	q.n[VW] = S; 
	q.n[VX] = Vx; 
	q.n[VY] = Vy; 
	q.n[VZ] = Vz;	
	return q;
}

//////////////////////////////////////////////////////////////////////////
// Convert Quaternion to rotation matrix
// Input is the Quaternion for conversion
void mat3::FromQuaternion(const Quaternion& q)
{
    
	double S = q.n[VW]; 
	double Vx = q.n[VX];
	double Vy = q.n[VY];
	double Vz = q.n[VZ];
	mat3 RotMat;
	RotMat[0][0] = 1 - 2*Vy*Vy - 2*Vz*Vz;
	RotMat[0][1] = 2*Vx*Vy - 2*S*Vz; 
	RotMat[0][2] = 2*Vx*Vz + 2*S*Vy; 
	RotMat[1][0] = 2*Vx*Vy + 2*S*Vz; 
	RotMat[1][1] = 1 - 2*Vx*Vx - 2*Vz*Vz; 
	RotMat[1][2] = 2*Vy*Vz - 2*S*Vx; 
	RotMat[2][0] = 2*Vx*Vz - 2*S*Vy; 
    RotMat[2][1] = 2*Vy*Vz+2*S*Vx; 
	RotMat[2][2] = 1-2*Vx*Vx - 2*Vy*Vy; 

	(*this) = RotMat; 
}

//========================================================================
// Quaternion

//////////////////////////////////////////////////////////////////////////
// Convert Quaternion to rotation matrix
// Return the corresponding matrix
mat3 Quaternion::ToRotation () const
{
	double S = n[VW]; 
	double Vx = n[VX];
	double Vy = n[VY];
	double Vz = n[VZ];

	mat3 RotMat;
	RotMat[0][0] = 1 - 2*Vy*Vy - 2*Vz*Vz;
	RotMat[0][1] = 2*Vx*Vy - 2*S*Vz; 
	RotMat[0][2] = 2*Vx*Vz + 2*S*Vy; 
	RotMat[1][0] = 2*Vx*Vy + 2*S*Vz; 
	RotMat[1][1] = 1 - 2*Vx*Vx - 2*Vz*Vz; 
	RotMat[1][2] = 2*Vy*Vz - 2*S*Vx; 
	RotMat[2][0] = 2*Vx*Vz - 2*S*Vy; 
    RotMat[2][1] = 2*Vy*Vz+2*S*Vx; 
	RotMat[2][2] = 1-2*Vx*Vx - 2*Vy*Vy; 

	return RotMat; 
}

//////////////////////////////////////////////////////////////////////////
// Convert a rotation matrix to Quaternion
// Input is the rotation matrix
// Return the corresponding Quaternion
// You need to handle the case when S is close to zero
void Quaternion::FromRotation (const mat3& rot)
{
	float tr=rot[0][0]+rot[1][1]+rot[2][2];
	float s;
	if(tr>0.0)
	{
		s=sqrt(1+rot[0][0]+rot[1][1]+rot[2][2]);
		n[VW]=s*0.5;
		s=0.5/s;
		
		n[VX]=(rot[2][1]-rot[1][2])*s;
		n[VY]=(rot[0][2]-rot[2][0])*s;
		n[VZ]=(rot[1][0]-rot[0][1])*s;
	}
	else
	{
		int order[3]={VY,VZ,VX};
		int i=VX;
		if(rot[VY][VY]>rot[VX][VX]) 
			i=VY;
		if(rot[VZ][VZ]>rot[i][i])
			i=VZ;
		int j=order[i];
		int k=order[j];

		s=sqrt(rot[i][i]-(rot[j][j]+rot[k][k])+1.0);
		n[i]=s*0.5;
		s=0.5/s;
		n[VW]=(rot[k][j]-rot[j][k])*s;
		n[j]=(rot[i][j]+rot[j][i])*s;
		n[k]=(rot[i][k]+rot[k][i])*s;

	}

}

//////////////////////////////////////////////////////////////////////////
// Slerp interpolation
// Input are: time = t, from Quaternion q0, to Quaternion q1
// Return the Quaternion that is the SLERP(t, q0, q1)
Quaternion Quaternion::Slerp(float t, const Quaternion& q0, const Quaternion& q1)
{
	//how can you make sure that t is between 0 and 1? 


	// Replace the following code with your code
	double CosOmiga = q0.n[VW]*q1.n[VW] + q0.n[VX]*q1.n[VX] + q0.n[VY]*q1.n[VY]+q0.n[VZ]*q1.n[VZ];
	double Omiga = acos(CosOmiga); 
	Quaternion q;
	q = sin((1-t)*Omiga)/sin(Omiga)*q0 + sin(Omiga*t)/sin(Omiga)*q1; 
	return q;
}

//////////////////////////////////////////////////////////////////////////
// Compute the intermediate Quaternion for SQUAD
// Input are three quaternions from frame(i-1), frame (i) and frame(i+1)
// Return the intermediate Quaternion for frame(i)
Quaternion Quaternion::Intermediate (const Quaternion& q0, const Quaternion& q1, const Quaternion& q2)
{
	// Replace the following code with your code
	Quaternion q;
	q = q1*q.Exp(-(q.Log(q1.Inverse()*q2) + q.Log(q1.Inverse()*q0))/4.0);      //exp(-(log(q1.Inverse()*q2)+log
	return q;
}

//////////////////////////////////////////////////////////////////////////
// Squad interpolation
// Compute the Squad interpolation
// Inputs are: time = t, from Quaternion q0, its intermediate Quaternion a,
//                       to Quaternion q1, its intermediate Quaternion b
// Return the Quaternion that is Squad(t, q0, a, b, q1)
Quaternion Quaternion::Squad(float t, const Quaternion& q0, const Quaternion& a, const Quaternion& b, const Quaternion& q1)
{
	// Replace the following code with your code
	Quaternion q;
	q = Slerp(2*t*(1-t), Slerp(t, q0, q1), Slerp(t, a, b)); 
	return q;
}

//=======================================================================
// Transform
//////////////////////////////////////////////////////////////////////////
// Inverse of a transformation
// Return the inverse of current transformation
Transform Transform::Inverse() const
{
	// Replace the following code with your code

	//think more about this part
	Transform tmp;
	mat3 rot; 
	rot = tmp.m_rotation; 
	rot = rot.Transpose(); 
	tmp.m_rotation = rot; 

    vec3 trans;
	trans = tmp.m_translation; 
	rot = -1 * rot; 
	trans = rot * trans; 
	tmp.m_translation = trans; 

	return tmp; 
}

//////////////////////////////////////////////////////////////////////////
// * operator for Transform
// Return the Transform that equals to transform1 * transform2
Transform operator * (const Transform& t1, const Transform& t2)
{
	// Replace the following code with your code
	Transform tmp;
	tmp.m_rotation = t1.m_rotation*t2.m_rotation; 
	tmp.m_translation = t1.m_rotation*t2.m_translation+t1.m_translation;  
	return tmp;
}
