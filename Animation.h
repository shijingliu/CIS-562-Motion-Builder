//////////////////////////////////////////////////////////////////////////
// Modified by
// Tiantian Liu
// 11/11/2011
// University of Pennsylvania


//////////////////////////////////////////////////////////////////////////
// Animation.h -- Header file for useful Classes for character animation
//
// Liming Zhao
// 11/02/2007
// University of Pennsylvania

#pragma  once

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

#include "Transformation.h"

using namespace std;

#ifndef ROTATIONORDER
#define ROTATIONORDER 
enum {ZXY, XYZ, ZYX, XZY, YXZ, YZX};
#endif

#ifndef EPSILONDIST
#define EPSILONDIST 0.1f
#endif

#ifndef MAX_ITER
#define MAX_ITER 5
#endif

class Joint;
class Skeleton;
class Frame;
class Motion;
class Player;

//////////////////////////////////////////////////////////////////////////
//
// Joint
//
//////////////////////////////////////////////////////////////////////////

class Joint{
public:
	Joint(void);
	Joint(const string& name);
	Joint(const int& id, const string& name);  
	~Joint(void);

	Joint* GetParent(void);
	unsigned int GetChildCount(void) const;
	Joint* GetChildAt(unsigned int index);

	void UpdateTransformation(bool bRecursive = false);

	void SetName(const string& name);
	void SetID(const int& id);
	void SetChannelCount(const unsigned int& count);
	void SetLocalTransform(const Transform& transform);
	void SetLocalTranslation(const vec3& translation);
	void SetLocalRotation(const mat3& rotation);
	void SetData(void* pData);

	const string& GetName(void) const;
	const int& GetID(void) const;
	const unsigned int& GetChannelCount(void) const;
	const Transform& GetLocalTransform(void) const;
	const vec3& GetLocalTranslation(void) const;
	const mat3& GetLocalRotation(void) const;
	const Transform& GetGlobalTransform(void) const;
	const vec3& GetGlobalTranslation(void) const;
	const mat3& GetGlobalRotation(void) const;
	void* GetData(void);

	static void AttachJoints(Joint* pParent, Joint* pChild);
	static void DetachJoints(Joint* pParent, Joint* pChild);

private:
	//Joint name
	string m_name;
	//Joint index number
	int m_id;
	//Channel count
	unsigned int m_channelCount;
	//Pointer to the parent joint
	Joint* m_pParent;
	//Vector of children joints
	vector<Joint*> m_children;
	//Transformation information
	Transform m_local;
	Transform m_global;
	//External Data
	void* m_pData;
};

//////////////////////////////////////////////////////////////////////////
//
// Skeleton
//
//////////////////////////////////////////////////////////////////////////

class Skeleton{
public:
	Skeleton(void);
	~Skeleton(void);

	bool LoadFromFile(ifstream& inFile);
	void SaveToFile(ofstream& outFile);
	Joint* GetJointByName(const string& name);
	Joint* GetJointByID(const unsigned int& id);
	Joint* GetRootJoint(void);
	const unsigned int& GetJointCount(void) const { return m_jointCount; }
	const int& GetSelectedJoint(void) { return m_selectedJoint; }	

	void UpdateFK(Joint* pRoot = NULL);
	void ReadFromFrame(Frame* pFrame);
	void WriteToFrame(Frame* pFrame);

	const vec3& GetGoalPosition(void) const;
	void SetGoalPosition(const vec3& pos);
	void SetSelectedJoint(const int& selectedJoint);	
	void SolveIK(void);

        friend class Player;

private:
	bool LoadFromFileRec(ifstream& inFile, Joint* pParent, string prefix);
	void SaveToFileRec(ofstream& outFile, Joint* pJoint, unsigned int level);
	void ReadFromFrameRec(Joint* pJoint, mat3* pRotationData);
	void WriteToFrameRec(Joint* pJoint, Frame* pFrame);
	vector<Joint*> m_joints;
	vector<Joint*> m_ikChain;
	unsigned int m_jointCount;
	Joint* m_pRoot;
	int m_selectedJoint;
	vec3 m_goalPosition;
};

//////////////////////////////////////////////////////////////////////////
//
// Frame
//
//////////////////////////////////////////////////////////////////////////

class Frame
{
public:
	Frame(void);
	~Frame(void);

	void SetJointCount(unsigned int jointCount);
	const unsigned int& GetJointCount(void) const;
	void LoadFromFile(ifstream& inFile, Skeleton* pSkeleton);	// Read from BVH file and assume ZXY rotation order
	void SaveToFile(ofstream& outFile, Skeleton* pSkeleton);	// Write to BVH file and assume ZXY rotation order
	void Clone(const Frame& frameData);		// Clone from given frame

	void SetRootTranslation(const vec3& translation);
	void SetJointRotation(const unsigned int& index, const mat3& rotation);
	void SetJointRotation(const unsigned int& index, const Quaternion& rotation);
	const vec3& GetRootTranslation(void) { return m_rootTranslation; }
	const mat3& GetJointRotation(const unsigned int& index) const { return m_rotationData[index]; }
	//const Quaternion& GetJointRotation(const unsigned int& index) {return m_quaternionData[index]; }

	static void Cubic(const Frame& frame0, const Frame& frame1, const Frame& frame2, const Frame& frame3, Frame& targetFrame, float fPerc);
	static void Slerp(const Frame& frame0, const Frame& frame1, Frame& targetFrame, float fPerc);
	static void Squad(const Frame& frame0, const Frame& s0,const Frame& s1, const Frame& frame1,
					  Frame& targetFrame, float fPerc);
	static void IntermediateFrame(const Frame& frame0, const Frame& frame1, const Frame& frame2, Frame& targetFrame);

        friend class Skeleton;

private:
	unsigned int m_jointCount;
	vec3 m_rootTranslation;
	mat3* m_rotationData;
	Quaternion* m_quaternionData;
};

//////////////////////////////////////////////////////////////////////////
//
// Motion
//
//////////////////////////////////////////////////////////////////////////
class Motion
{
public:
	Motion(void);
	~Motion(void);
	
	bool LoadFromFile(ifstream& inFile, Skeleton* pSkeleton);	// Read from BVH file
	void SaveToFile(ofstream& outFile, Skeleton* pSkeleton);	// Write to BVH file	
	void Clone(const Motion& motionData);	// Clone the motion

	Frame* GetFrameAt(const unsigned int& index);							// Get frame at index, allowing reading and writing
	Frame* GetCurrentFrame(void);
	void SetCurrentFrameCount(const unsigned int& currentFrame) { m_currentFrame = currentFrame; }
	const unsigned int& GetCurrentFrameCount(void) { return m_currentFrame; }
	void Update(int deltaFrame);
	const unsigned int& GetJointCount(void) const {return m_jointCount;}	// Get joint count
	const unsigned int& GetFrameCount(void) const {return m_frameCount;}	// Get frame count
	const string& GetName(void) const {return m_name;}	// Get name
	void SetName(const string& name) {m_name = name;}	// Set name	

	// Create an animation and allocate spaces
	static Motion* CreateMotion(unsigned int jointCount, unsigned int frameCount, string name);

private:
	void FreeSpace(void);
	void AllocateSpace(void);

	vector<Frame*> m_keyFrames;
	unsigned int m_jointCount;
	unsigned int m_frameCount;
	unsigned int m_currentFrame;
	string m_name;

        friend class Player;
};

//////////////////////////////////////////////////////////////////////////
//
// Player
//
//////////////////////////////////////////////////////////////////////////
enum{eMotion1 = 1, eMotion2, eMotion3};
enum{eCUBIC, eSLERP, eSQUAD};

class Player
{
public:
	Player(void);
	~Player(void);

	bool LoadBVHFile(ifstream& inFile);
	bool SaveBVHFile(ofstream& outFile);
	void UpdateFrame(void);
	void SaveFrame(void);

	bool IsValid(void);
	bool BlendingReady(unsigned int& frameCount0, unsigned int& frameCount1);
	void SetCurrentMotion(const unsigned int& index);
	const unsigned int& GetCurrentMotionIndex(void) const { return m_motionIndex; }
	Motion* GetCurrentMotion(void);
	Skeleton* GetSkeleton(void);

	void Blend(const unsigned int& startFrame, const unsigned int& endFrame, const unsigned int& interpFrame, const unsigned int& type);
	
	const vec3& GetGoalPosition(void) const;
	void SetGoalPosition(const vec3& pos);
	const unsigned int& GetMode(void) const;
	void SetMode(const unsigned int& mode);
	void SolveIK(void);

private:
	Skeleton m_skeleton;
	Motion m_motion1, m_motion2, m_motion3;
	Motion *m_currentMotion;
	unsigned int m_motionIndex;
};
