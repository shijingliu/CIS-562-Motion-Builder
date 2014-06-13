//////////////////////////////////////////////////////////////////////////
// Modified by
// Tiantian Liu
// 11/11/2011
// University of Pennsylvania

//////////////////////////////////////////////////////////////////////////
// Animation.cpp -- Source file for useful Classes for character animation
//
// Liming Zhao
// 11/02/2007
// University of Pennsylvania

#include "Animation.h"

//////////////////////////////////////////////////////////////////////////
//
// Joint
//
//////////////////////////////////////////////////////////////////////////



// Constructors
Joint::Joint(void)
{
	m_name = "";
	m_id = -1;
	m_pParent = NULL;
	m_pData = NULL;
}

Joint::Joint(const int& id, const string& name)
{
	m_id = id;
	m_name = name;
	m_pParent = NULL;
	m_pData = NULL;
}

Joint::Joint(const string& name)
{
	m_name = name;
	m_pParent = NULL;
	m_pData = NULL;
}

// Destructor
Joint::~Joint(void)
{
}

// Member functions  
Joint* Joint::GetParent(void)
{
	return m_pParent;
}

unsigned int Joint::GetChildCount(void) const
{
	return unsigned int (m_children.size());
}

Joint* Joint::GetChildAt(unsigned int index)
{
	return m_children[index];
}

void Joint::SetName(const string& name)
{
	m_name = name;
}

void Joint::SetID(const int& id)
{
	m_id = id;
}

void Joint::SetChannelCount(const unsigned int& count)
{
	m_channelCount = count;
}

void Joint::SetLocalTransform(const Transform& transform)
{
	m_local = transform;
}

void Joint::SetLocalTranslation(const vec3& translation)
{
	m_local.m_translation = translation;
}

void Joint::SetLocalRotation(const mat3& rotation)
{
	m_local.m_rotation = rotation;
}

void Joint::SetData(void* pData)
{
	m_pData = pData;
}

const string& Joint::GetName(void) const
{
	return m_name;
}

const int& Joint::GetID(void) const
{
	return m_id;
}

const unsigned int& Joint::GetChannelCount(void) const
{
	return m_channelCount;
}

const Transform& Joint::GetLocalTransform(void) const
{
	return m_local;
}

const vec3& Joint::GetLocalTranslation(void) const
{
	return m_local.m_translation;
}

const mat3& Joint::GetLocalRotation(void) const
{
	return m_local.m_rotation;
}

const Transform& Joint::GetGlobalTransform(void) const
{
	return m_global;
}

const vec3& Joint::GetGlobalTranslation(void) const
{
	return m_global.m_translation;
}

const mat3& Joint::GetGlobalRotation(void) const
{
	return m_global.m_rotation;
}

void* Joint::GetData(void)
{
	return m_pData;
}

void Joint::AttachJoints(Joint* pParent, Joint* pChild)
{
	if (pChild)
	{
		Joint* pOldParent = pChild->m_pParent;
		if (pOldParent)
		{
			// erase the child from old parent's children list
			vector<Joint*>::iterator iter;
			for (iter = pOldParent->m_children.begin(); iter != pOldParent->m_children.end(); iter++)
			{
				if ((*iter) == pChild)
				{
					iter = pOldParent->m_children.erase(iter);
				}
			}
		}
		// Set the new parent
		pChild->m_pParent = pParent;
		// Add child to new parent's children list
		if (pParent)
		{
			pParent->m_children.push_back(pChild);
		}
	}
}

void Joint::DetachJoints(Joint* pParent, Joint* pChild)
{
	if (pChild && pChild->m_pParent == pParent)
	{
		if (pParent)
		{
			// erase the child from parent's children list
			vector<Joint*>::iterator iter;
			for (iter = pParent->m_children.begin(); iter != pParent->m_children.end(); iter++)
			{
				if ((*iter) == pChild)
				{
					iter = pParent->m_children.erase(iter);
				}
			}
		}
		pChild->m_pParent = NULL;
	}
}

//////////////////////////////////////////////////////////////////////////
//
// Skeleton
//
//////////////////////////////////////////////////////////////////////////

Skeleton::Skeleton(void)
{
	m_pRoot = NULL;
	m_jointCount = 0;
	m_selectedJoint = -1;
}

Skeleton::~Skeleton(void)
{	
	m_ikChain.clear();
	for(vector<Joint*>::iterator iter = m_joints.begin(); iter != m_joints.end(); iter++)
	{
		delete (*iter);
	}
	m_joints.clear();
	
}

bool Skeleton::LoadFromFile(ifstream& inFile)
{
	if (m_joints.size() > 0)
	{
		m_ikChain.clear();
		for(vector<Joint*>::iterator iter = m_joints.begin(); iter != m_joints.end(); iter++)
		{
			delete (*iter);
		}
		m_joints.clear();
		m_jointCount = 0;
		m_pRoot = NULL;
		m_selectedJoint = -1;
		m_goalPosition = vec3Zero;
	}

	string readString, jointname;
	unsigned int channelCount;
		
	inFile >> readString;
	if (readString != "HIERARCHY")
		return false;
	inFile >> readString;
	if (readString != "ROOT")
		return false;
	inFile.get(); //" "
	getline(inFile, jointname);// joint name
	Joint* joint = new Joint(jointname);
	m_joints.push_back(joint);
	joint->SetID(unsigned(m_joints.size() - 1));
	joint->SetChannelCount(6);
	m_pRoot = joint;
	inFile >> readString; // "{"
	inFile >> readString; // "OFFSET"
	getline(inFile, readString);	// Root translation values
	inFile >> readString;
	if (readString != "CHANNELS")
		return false;
	inFile >> channelCount;
	joint->SetChannelCount(channelCount);
	getline(inFile, readString);	// " Xposition Yposition Zposition Zrotation Xrotation Yrotation"
	inFile >> readString;
	while(readString != "}") 
	{		
		if (!LoadFromFileRec(inFile, joint, readString))
		{
			return false;
		}
		inFile >> readString;
	}
	if (readString != "}")
		return false;

	m_jointCount = unsigned int(m_joints.size());
	return true;
}

bool Skeleton::LoadFromFileRec(ifstream &inFile, Joint *pParent, string prefix)
{
	string readString, jointname;
	vec3 offsets;
	unsigned int channelCount;
	if (prefix == "JOINT")
	{
		inFile.get(); //" "
		getline(inFile, jointname);// joint name
		Joint* joint = new Joint(jointname);
		m_joints.push_back(joint);
		joint->SetID(unsigned(m_joints.size() - 1));	
		Joint::AttachJoints(pParent, joint);
		inFile >> readString; // "{"
		inFile >> readString; // "OFFSET"
		inFile >> offsets[0] >> offsets[1] >> offsets[2];
		joint->SetLocalTranslation(offsets);
		inFile >> readString; // "CHANNELS"
		inFile >> channelCount;
		joint->SetChannelCount(channelCount);
		getline(inFile, readString);// " Zrotation Xrotation Yrotation"
		inFile >> readString; // "Joint" or "}" or "End"
		while (readString != "}")
		{
			if (LoadFromFileRec(inFile, joint, readString) == false)
				return false;
			inFile >> readString; // "Joint" or "}" or "End"
		}
		return true;
	}else if (prefix == "End")
	{	
		inFile.get(); //" "
		getline(inFile, jointname);// joint name
		Joint* joint = new Joint(jointname);
		m_joints.push_back(joint);
		joint->SetID(unsigned(m_joints.size() - 1));
		joint->SetChannelCount(0);
		Joint::AttachJoints(pParent, joint);
		inFile >> readString; // "{"
		inFile >> readString; // "OFFSET"
		inFile >> offsets[0] >> offsets[1] >> offsets[2];
		joint->SetLocalTranslation(offsets);
		inFile >> readString; // "}"
		return true;
	}else
		return false;
}

void Skeleton::SaveToFile(ofstream& outFile)
{
	outFile << "HIERARCHY" << endl;
	outFile << "ROOT " << m_pRoot->GetName() << endl;
	outFile << "{" << endl;
	outFile << "\tOFFSET 0.00 0.00 0.00" << endl;
	outFile << "\tCHANNELS " << m_pRoot->GetChannelCount() << " Xposition Yposition Zposition Zrotation Xrotation Yrotation" << endl;
	unsigned int childCount = m_pRoot->GetChildCount();
	for (unsigned int i = 0; i < childCount; i++)
	{
		Joint* pChild = m_pRoot->GetChildAt(i);
		SaveToFileRec(outFile, pChild, 1);
	}
	outFile << "}" << endl;
}

void Skeleton::SaveToFileRec(ofstream& outFile, Joint* pJoint, unsigned int level)
{
	string indentation = "";
	for (unsigned int i = 0; i < level; i++)
	{
		indentation += "\t";
	}

	if (pJoint->GetChannelCount() == 3)
	{
		outFile << indentation << "JOINT " << pJoint->GetName() << endl;
		outFile << indentation << "{" << endl;
		const vec3& offsets = pJoint->GetLocalTranslation();
		outFile << indentation << "\tOFFSET " << offsets[0] << " " << offsets[1] << " " << offsets[2] << endl;
		outFile << indentation << "\tCHANNELS " << pJoint->GetChannelCount() << " Zrotation Xrotation Yrotation" << endl;
	}else
	{
		outFile << indentation << "End " << pJoint->GetName() << endl;
		outFile << indentation << "{" << endl;
		const vec3& offsets = pJoint->GetLocalTranslation();
		outFile << indentation << "\tOFFSET " << offsets[0] << " " << offsets[1] << " " << offsets[2] << endl;		
	}
	unsigned int childCount = pJoint->GetChildCount();
	for (unsigned int i = 0; i < childCount; i++)
	{
		Joint* pChild = pJoint->GetChildAt(i);
		SaveToFileRec(outFile, pChild, level + 1);
	}
	outFile << indentation << "}" << endl;
}

Joint* Skeleton::GetJointByName(const string& name)
{
	vector<Joint*>::const_iterator iter;
	for (iter = m_joints.begin(); iter != m_joints.end(); iter++)
	{
		if (name == ((*iter)->GetName()))
			return (*iter);
	}
	return NULL;
}

Joint* Skeleton::GetJointByID(const unsigned int& id)
{
	if (id >= 0 && id < m_jointCount)
	{
		return m_joints[id];
	}else
		return NULL;
}

Joint* Skeleton::GetRootJoint(void)
{
	return m_pRoot;
}

void Skeleton::ReadFromFrame(Frame* pFrame)
{
	if (m_jointCount != pFrame->GetJointCount())
		return;
	m_pRoot->SetLocalTranslation(pFrame->m_rootTranslation);
	m_pRoot->SetLocalRotation(pFrame->m_rotationData[m_pRoot->GetID()]);
	unsigned int totalChildren = m_pRoot->GetChildCount();
	for (unsigned int i = 0; i < totalChildren; i++)
	{
		ReadFromFrameRec(m_pRoot->GetChildAt(i), pFrame->m_rotationData);
	}
	UpdateFK(m_pRoot);
}

void Skeleton::ReadFromFrameRec(Joint* pJoint, mat3* pRotationData)
{
	pJoint->SetLocalRotation(pRotationData[pJoint->GetID()]);
	unsigned int totalChildren = pJoint->GetChildCount();
	for (unsigned int i = 0; i < totalChildren; i++)
	{
		ReadFromFrameRec(pJoint->GetChildAt(i), pRotationData);
	}
}

void Skeleton::WriteToFrame(Frame* pFrame)
{
	if (m_jointCount != pFrame->GetJointCount())
		return;
	pFrame->m_rootTranslation = m_pRoot->GetLocalTranslation();
	pFrame->SetJointRotation(m_pRoot->GetID(), m_pRoot->GetLocalRotation());
	unsigned int totalChildren = m_pRoot->GetChildCount();
	for (unsigned int i = 0; i < totalChildren; i++)
	{
		WriteToFrameRec(m_pRoot->GetChildAt(i), pFrame);
	}
}

void Skeleton::WriteToFrameRec(Joint* pJoint, Frame* pFrame)
{
	pFrame->SetJointRotation(pJoint->GetID(), pJoint->GetLocalRotation());
	unsigned int totalChildren = pJoint->GetChildCount();
	for (unsigned int i = 0; i < totalChildren; i++)
	{
		WriteToFrameRec(pJoint->GetChildAt(i), pFrame);
	}
}

const vec3& Skeleton::GetGoalPosition(void) const
{
	return m_goalPosition;
}

void Skeleton::SetGoalPosition(const vec3& pos)
{
	m_goalPosition = pos;
}

void Skeleton::SetSelectedJoint(const int& selectedJoint)
{
	m_selectedJoint = selectedJoint;
	m_ikChain.clear();
	if (m_selectedJoint >= 0 && m_selectedJoint < int(m_jointCount))
	{		
		Joint* tmp = m_joints[selectedJoint];
		m_goalPosition = tmp->GetGlobalTranslation();
		if (tmp != NULL && tmp->GetParent() != NULL)
		{
			// terminate at the joint before root joint, so that root will not change during IK	
			tmp = tmp->GetParent();
			while(tmp->GetParent() != NULL)
			{
				m_ikChain.push_back(tmp);
				tmp = tmp->GetParent();
			}
		}
	}
	
}


//////////////////////////////////////////////////////////////////////////
//
// Frame
//
//////////////////////////////////////////////////////////////////////////
Frame::Frame(void)
{
	m_jointCount = 0;
	m_rotationData = NULL;
	m_quaternionData = NULL;
}

Frame::~Frame(void)
{
	if (m_rotationData)
		delete[] m_rotationData;
	if (m_quaternionData)
	{
		delete[] m_quaternionData;
	}
}

void Frame::SetJointCount(unsigned int jointCount)
{
	if (m_jointCount != jointCount)
	{
		m_jointCount = jointCount;
		if (m_rotationData)
		{
			delete[] m_rotationData;
		}
		m_rotationData = new mat3[jointCount];
		if (m_quaternionData)
		{
			delete[] m_quaternionData;
		}
		m_quaternionData = new Quaternion[jointCount];
	}
}

const unsigned int& Frame::GetJointCount(void) const
{
	return m_jointCount;
}

void Frame::LoadFromFile(ifstream& inFile, Skeleton* pSkeleton)
{
	float tx, ty, tz, ry, rx, rz;
	inFile >> tx >> ty >> tz;
	m_rootTranslation = vec3(tx, ty, tz);
	for (unsigned int i = 0; i < m_jointCount; i++)
	{
		Joint* pJoint = pSkeleton->GetJointByID(i);
		if (pJoint->GetChannelCount() > 0)
		{
			inFile >> rz >> rx >> ry;
		}else
		{
			rz = rx = ry = 0.0f;
		}
		m_rotationData[i].FromEulerAnglesZXY(vec3(rx, ry, rz) * Deg2Rad);
		m_quaternionData[i].FromRotation(m_rotationData[i]);
	}
}

void Frame::SaveToFile(ofstream& outFile, Skeleton* pSkeleton)
{
	outFile << setprecision(6);
	outFile << m_rootTranslation[0] << "\t" << m_rootTranslation[1] << "\t" << m_rootTranslation[2];
	vec3 angles;
	for (unsigned int i = 0; i < m_jointCount; i++)
	{
		Joint* pJoint = pSkeleton->GetJointByID(i);
		if (pJoint->GetChannelCount() > 0)
		{
			m_rotationData[i].ToEulerAnglesZXY(angles);
			angles *= Rad2Deg;
			outFile << "\t" << angles[VZ] << "\t" << angles[VX] << "\t" << angles[VY];
		}
	}
	outFile << endl;
}

void Frame::Clone(const Frame& frameData)
{
	SetJointCount(frameData.m_jointCount);
	m_rootTranslation = frameData.m_rootTranslation;
	for (unsigned int i = 0; i < m_jointCount; i++)
	{
		m_rotationData[i] = frameData.m_rotationData[i];
		m_quaternionData[i] = frameData.m_quaternionData[i];
	}
}

void Frame::SetRootTranslation(const vec3& translation)
{
	m_rootTranslation = translation;
}

void Frame::SetJointRotation(const unsigned int& index, const mat3& rotation)
{
	assert(index >= 0 && index < m_jointCount);
	m_rotationData[index] = rotation;
	m_quaternionData[index].FromRotation(rotation);
}

void Frame::SetJointRotation(const unsigned int& index, const Quaternion& rotation)
{
	assert(index >= 0 && index < m_jointCount);	
	m_quaternionData[index] = rotation;
	m_rotationData[index].FromQuaternion(rotation);
}
//////////////////////////////////////////////////////////////////////////
//
vec3 CubicVec3(vec3 &d1, vec3 &d2, vec3 &d3, vec3 &d4, float t)
{
	vec3 a = d2; 
	vec3 b = d2 - d1; 
	vec3 c = (d3 - d2) * 3 - (d2 - d1) * 2 - (d4 - d3); 
	vec3 d = (d2 - d3) * 2 + d2 - d1 + d4 - d3;
	return a + b * t + c * t * t + d * t * t * t;
}


//////////////////////////////////////////////////////////////////////////
//
// Motion
//
//////////////////////////////////////////////////////////////////////////

Motion::Motion(void)
{
	m_jointCount = 0;
	m_name = "unknown";
	m_frameCount = 0;
	m_currentFrame = 0;
}

Motion::~Motion(void)
{
	FreeSpace();
}

bool Motion::LoadFromFile(ifstream& inFile, Skeleton* pSkeleton)
{
	string readString;
	unsigned int frameCount;
	inFile >> readString;
	if (readString != "MOTION")
		return false;
	inFile >> readString;
	if(readString != "Frames:") 
		return false;
	inFile >> frameCount;
	inFile >> readString; // "Frame"
	getline(inFile, readString); // " Time: 0.033333"

	FreeSpace();
	m_jointCount = pSkeleton->GetJointCount();
	m_frameCount = frameCount;
	AllocateSpace();
	
	for (vector<Frame*>::const_iterator iter = m_keyFrames.begin(); iter != m_keyFrames.end(); iter++)
	{
		(*iter)->LoadFromFile(inFile, pSkeleton);
	}
	m_currentFrame = 0;
	return true;
}  

void Motion::SaveToFile(ofstream& outFile, Skeleton* pSkeleton)
{
	outFile << "MOTION" << endl;
	outFile << "Frames: " << m_frameCount << endl;
	outFile << "Frame Time: 0.033333" << endl;
	for (vector<Frame*>::const_iterator iter = m_keyFrames.begin(); iter != m_keyFrames.end(); iter++)
	{
		(*iter)->SaveToFile(outFile, pSkeleton);
	}
}

void Motion::FreeSpace(void)
{
	for (vector<Frame*>::const_iterator iter = m_keyFrames.begin(); iter != m_keyFrames.end(); iter++)
	{
		delete (*iter);
	}
	m_keyFrames.clear();
	m_frameCount = 0;
	m_jointCount = 0;
}

void Motion::AllocateSpace(void)
{
	Frame* pFrameData;
	for (unsigned int i = 0; i < m_frameCount; i++)
	{
		pFrameData = new Frame();
		pFrameData->SetJointCount(m_jointCount);
		m_keyFrames.push_back(pFrameData);
	}
}


void Motion::Clone(const Motion& motionData)
{
	// Free old frames
	if (m_keyFrames.size() > 0)
	{
		FreeSpace();
	}
	
	// Clone info
	m_name = motionData.m_name;
	m_jointCount = motionData.m_jointCount;
	m_frameCount = motionData.m_frameCount;

	//Clone animation data
	Frame* pFrameData;
	for (unsigned int i = 0; i < m_frameCount; i++)
	{
		pFrameData = new Frame();
		pFrameData->Clone(*(motionData.m_keyFrames[i]));
		m_keyFrames.push_back(pFrameData);
	}
}

Frame* Motion::GetFrameAt(const unsigned int& index)
{
	assert(index >= 0 && index < m_frameCount);
	return m_keyFrames[index];
}

Frame* Motion::GetCurrentFrame(void)
{
	return m_keyFrames[m_currentFrame];
}

Motion* Motion::CreateMotion(unsigned int jointCount, unsigned int frameCount, string name)
{
	Motion* pMotion = new Motion();
	pMotion->m_jointCount = jointCount;
	pMotion->m_frameCount = frameCount;
	pMotion->m_name = name;
	pMotion->AllocateSpace();
	return pMotion;
}

void Motion::Update( int deltaFrame )
{
	int num = int(m_currentFrame) + deltaFrame;
	while (num >= int(m_frameCount))
		num -= m_frameCount;
	while (num < 0)
		num += m_frameCount;
	m_currentFrame = unsigned int(num);
}
//////////////////////////////////////////////////////////////////////////
//
// Player
//
//////////////////////////////////////////////////////////////////////////

Player::Player(void)
{
	m_currentMotion = &m_motion1;
	m_motionIndex = eMotion1;
}

Player::~Player(void)
{
}

bool Player::LoadBVHFile(ifstream& inFile)
{
	if (m_skeleton.LoadFromFile(inFile) == false)
		return false;
	if (m_currentMotion->LoadFromFile(inFile, &m_skeleton) == false)
		return false;
	return true;
}

bool Player::SaveBVHFile(ofstream& outFile)
{
	m_skeleton.SaveToFile(outFile);
	m_currentMotion->SaveToFile(outFile, &m_skeleton);
	return true;
}

void Player::UpdateFrame(void)
{
	if (IsValid())
	{
		m_skeleton.ReadFromFrame(m_currentMotion->GetCurrentFrame());
	}
}

void Player::SetCurrentMotion(const unsigned int& index)
{
	m_motionIndex = index;
	switch(index)
	{
	case eMotion1:
		m_currentMotion = &m_motion1;
		break;
	case eMotion2:
		m_currentMotion = &m_motion2;
	    break;
	case eMotion3:
		m_currentMotion = &m_motion3;
	    break;
	}
}
Motion* Player::GetCurrentMotion(void)
{
	return m_currentMotion;
}

Skeleton* Player::GetSkeleton(void)
{
	return &m_skeleton;
}

bool Player::IsValid(void)
{
	return (m_skeleton.GetJointCount() > 0 && m_currentMotion && m_currentMotion->GetFrameCount() > 0);
}

void Player::SaveFrame(void)
{
	m_skeleton.WriteToFrame(m_currentMotion->GetCurrentFrame());
}

bool Player::BlendingReady(unsigned int& frameCount0, unsigned int& frameCount1)
{
	frameCount0 = m_motion1.GetFrameCount();
	frameCount1 = m_motion2.GetFrameCount();
	return (m_motion1.GetJointCount() == m_motion2.GetJointCount() && m_motion1.GetFrameCount() > 0 && m_motion2.GetFrameCount() >0);
}

const vec3& Player::GetGoalPosition(void) const
{
	return m_skeleton.m_goalPosition;
}

void Player::SetGoalPosition(const vec3& pos)
{
	m_skeleton.m_goalPosition = pos;
}

void Player::SolveIK(void)
{
	m_skeleton.SolveIK();
}


//************************************************************************
// Functions that you need to implement
//************************************************************************

//========================================================================
// Joint

//////////////////////////////////////////////////////////////////////////
// UpdateTransformation
// Input: whether to propagate the transformation update to its children joints
// If bRecursive == true, you need to update the all children joints transformations
// This basically implements the forward kinematics from current joint
void Joint::UpdateTransformation(bool bRecursive)
{
	if(this->GetID()!=0)
	{
		m_global=m_pParent->GetGlobalTransform()*m_local;
	}
	else
	{
		m_global=m_local;
	}
	if (bRecursive)
	{
		vector<Joint*>::const_iterator iter;
		for (iter = m_children.begin(); iter != m_children.end(); iter++)
		{
			(*iter)->UpdateTransformation(true);
		}
	}
}

//========================================================================
// Skeleton
//////////////////////////////////////////////////////////////////////////
// Update forward kinematics from given Joint pRoot
// If pRoot == NULL, update from root of the skeleton, which is m_pRoot 
void Skeleton::UpdateFK(Joint* pRoot)
{
	if(pRoot!=NULL)
	{
		pRoot->UpdateTransformation(true);
	}
	else
	{
		m_pRoot->UpdateTransformation(true);
	}
}

//////////////////////////////////////////////////////////////////////////
// Solve Inverse Kinematics using CCD approach
// When this functions is called, the following are given
// IK chain from selected Joint to root is store in m_ikChain. 
//		However, the first Joint is the parent of the selected Joint, the last Joint is the child of m_pRoot.
//		Basically you are operating the IK chain without moving the Skeleton root.
// End-effector(the selected Joint) is given as Joint* pEndEffector = m_joints[m_selectedJoint];
// Goal position is in m_goalPosition in world coordinates
// Maximum iteration is defined in MAX_ITER
// Please follow the structure of given code        
void Skeleton::SolveIK(void)
{
	if (m_ikChain.size() > 0) // If there is at east some Joints in the IK chain for manipulation
	{
		Joint* pEndEffector = m_joints[m_selectedJoint]; // Get the end-effector
		vec3 dist = m_goalPosition - pEndEffector->GetGlobalTranslation(); // Compute distance error
		//declare the last chain here 
		//Joint* last_chain=pEndEffector;
		if (dist.Length() > EPSILONDIST) // There is need for CCD
		{
			// You might want to do some preprocessing here
			for (unsigned int i = 0; i < MAX_ITER; i++)
			{
				Joint* pJoint;
				for (vector<Joint*>::const_iterator iter = m_ikChain.begin(); iter != m_ikChain.end(); ++iter)
				{
					// Let pJoint be the current Joint on the IK chain for handling      
					pJoint = *iter;
					
					vec3 r= (pJoint->GetGlobalRotation().Transpose())*(-pJoint->GetGlobalTranslation()+pEndEffector->GetGlobalTranslation());
					dist = (pJoint->GetGlobalRotation().Transpose())*dist; 
					
					vec3 axis=r.Cross(dist);
					axis = axis.Normalize(); 
					if(axis.Length() <0.001)
						continue; 

					float theta = 0.0; 
					float c = 0.1; 

					for(int i = 0; i< 10; ++i)
					{
						theta= (c*axis.Length()/(r.Length()*r.Length()+r*dist));
						mat3 T=mat3::Rotation3DRad(axis,theta); 

						pJoint->SetLocalRotation(T*pJoint->GetLocalRotation()); 	
						UpdateFK(pJoint);
					}
		     
					// Check error, if it is close enough, terminate the CCD iteration
					dist = m_goalPosition - pEndEffector->GetGlobalTranslation();
					if (dist.Length() < EPSILONDIST)
						return;
				}
			}
		}
	}
}

//========================================================================
// Frame
//////////////////////////////////////////////////////////////////////////
// Cubic interpolation from frame[start] to frame[end] using frame[start -1] and frame[end + 1] to help to compute the derivatives for Cubic interpolation
// Inputs are: frame0 = frame[start - 1], frame1 = frame[start], frame2 = frame[end], frame3 = frame[end + 1],
//		targetFrame holds the result, fPerc is between [0,1] for interpolation
// The root translation is interpolated as 3D vector cubic interpolation
// The joint rotations are interpolated using cubic interpolation on Euler angles. 
//		You need to convert the rotation matrix into Euler angles in ZXY order, and treat them as 3D vector and apply cubic interpolation.
//		Finally, map the interpolated Euler angles back to rotation matrix
// Remember that every frame carries both rotation matrix and quaternion,
//		In order to keep data consistency, you need to update the Quaternion when you modify rotation matrix and vice versa.
void Frame::Cubic(const Frame& frame0, const Frame& frame1, const Frame& frame2, const Frame& frame3, Frame& targetFrame, float fPerc)
{
	// Replace the following code with your code 
	vec3 angles;
	vec3 angle0;
	vec3 angle1; 
	vec3 angle2; 
	vec3 angle3; 

	vec3 pStart = frame1.m_rootTranslation; 
	vec3 pEnd = frame2.m_rootTranslation; 
	vec3 vStart = frame1.m_rootTranslation - frame0.m_rootTranslation; 
	vec3 vEnd = frame3.m_rootTranslation - frame2.m_rootTranslation; 
	targetFrame.m_rootTranslation = (1-3*fPerc*fPerc+2*fPerc*fPerc*fPerc)*pStart + fPerc*fPerc*(3-2*fPerc)*pEnd + fPerc*(fPerc-1)*(fPerc-1)*vStart+fPerc*fPerc*(fPerc-1)*vEnd; 

	for(unsigned int i = 0; i < frame0.m_jointCount; i++)
	{
		frame0.m_rotationData[i].ToEulerAnglesZXY(angle0);
		frame1.m_rotationData[i].ToEulerAnglesZXY(angle1);
		frame2.m_rotationData[i].ToEulerAnglesZXY(angle2);
		frame3.m_rotationData[i].ToEulerAnglesZXY(angle3);

	    vec3 vStAngle = angle1 - angle0; 
		vec3 vEnAngle = angle3 - angle2; 
		angles = (1-3*fPerc*fPerc+2*fPerc*fPerc*fPerc)*angle1 + fPerc*fPerc*(3-2*fPerc)*angle2 + fPerc*(fPerc-1)*(fPerc-1)*vStAngle+fPerc*fPerc*(fPerc-1)*vEnAngle; 
	    targetFrame.m_rotationData[i].FromEulerAnglesZXY(angles);
		targetFrame.m_quaternionData[i].FromRotation(targetFrame.m_rotationData[i]);
	}
}


//////////////////////////////////////////////////////////////////////////
// Slerp interpolation from frame[start] to frame[end]
// Inputs are: frame0 = frame[start], frame1 = frame[end], targetFrame holds the result, fPerc is between [0,1] for interpolation
// The root translation is interpolated as 3D vector linear interpolation
// The joint rotations are interpolated using Quaternion::Slerp()
// Remember that every frame carries both rotation matrix and quaternion,
//		In order to keep data consistency, you need to update the Quaternion when you modify rotation matrix and vice versa.
void Frame::Slerp(const Frame& frame0, const Frame& frame1, Frame& targetFrame, float fPerc)
{
	vec3 translate1 = frame0.m_rootTranslation; 
	vec3 translate2 = frame1.m_rootTranslation; 
	vec3 newTrans = (1-fPerc)*translate1 + fPerc*translate2; 
	targetFrame.m_rootTranslation = newTrans; 
	for(unsigned int i = 0; i< frame0.m_jointCount; i++)
	{
		targetFrame.m_quaternionData[i]=Quaternion::Slerp(fPerc,frame0.m_quaternionData[i],frame1.m_quaternionData[i]);
		targetFrame.m_rotationData[i]=targetFrame.m_quaternionData[i].ToRotation();
	}
}

//////////////////////////////////////////////////////////////////////////
// Compute the intermediate frame for Squad interpolation on frame[i]
// Inputs are: frame0 = frame[i - 1], frame1 = frame[i], frame2 = frame[i + 1], targetFrame holds the result
// Use Quaternion::Intermediate()
void Frame::IntermediateFrame(const Frame& frame0, const Frame& frame1, const Frame& frame2, Frame& targetFrame)
{
	// Add your code here
	Quaternion q0;
	Quaternion q1;
	Quaternion q2;
	for (unsigned int i = 0; i < frame0.m_jointCount; i++)
	{
		q0=frame0.m_rotationData[i].ToQuaternion();
		q1=frame1.m_rotationData[i].ToQuaternion();
		q2=frame2.m_rotationData[i].ToQuaternion();

		targetFrame.m_quaternionData[i]=Quaternion::Intermediate(q0,q1,q2);
		targetFrame.m_rotationData[i]=targetFrame.m_quaternionData[i].ToRotation();
	} 
}

//////////////////////////////////////////////////////////////////////////
// Squad interpolation from frame[start] to frame[end]
// Inputs are: frame0 = frame[start], s0  = intermediate frame for frame[start],
//			   s1 = intermediate frame for frame[end], frame1 = frame[end],
//			   targetFrame holds the result, fPerc is between [0,1] for interpolation
// The root translation is interpolated as 3D vector linear interpolation
// The joint rotations are interpolated using Quaternion::Squad()
// Remember that every frame carries both rotation matrix and quaternion,
//		In order to keep data consistency, you need to update the Quaternion when you modify rotation matrix and vice versa.
void Frame::Squad(const Frame& frame0, const Frame& s0,const Frame& s1, const Frame& frame1,
				  Frame& targetFrame, float fPerc)
{
	vec3 translate1 = frame0.m_rootTranslation; 
	vec3 translate2 = frame1.m_rootTranslation; 
	vec3 newTrans = (1-fPerc)*translate1 + fPerc*translate2; 
	targetFrame.m_rootTranslation = newTrans; 
	for (unsigned int i = 0; i < frame0.m_jointCount; i++)
	{
		// Slerp the joint rotations
		targetFrame.m_quaternionData[i]=Quaternion::Squad(fPerc,frame0.m_quaternionData[i],s0.m_quaternionData[i],s1.m_quaternionData[i],frame1.m_quaternionData[i]);
		targetFrame.m_rotationData[i]=targetFrame.m_quaternionData[i].ToRotation();
		// Make sure of the data consistency
	}
}


//========================================================================
// Player
//////////////////////////////////////////////////////////////////////////
// Blending from startFrame of m_motion1 into endFrame of m_motion2 with interpFrame number of interpolation frames and interpolation type is in type
// Follow the instructions below to complete this function
// Write some member functions if they help you to solve the problem    

void Player::Blend( const unsigned int& startFrame, const unsigned int& endFrame, const unsigned int& interpFrame, const unsigned int& type )
{
    // Compute total frame number for the target motion
    //unsigned int totalFrame = 1; // replace with your code
	unsigned int totalFrame = startFrame + interpFrame + (m_motion2.GetFrameCount() - endFrame) + 1;
	
    // Release and allocation space for result motion m_motion3
    m_motion3.FreeSpace();
    m_motion3.m_jointCount = m_skeleton.GetJointCount();
    m_motion3.m_frameCount = totalFrame;
    m_motion3.AllocateSpace();

    // Copy frames from m_motion1 to m_motion3
    // m_motion3[0, ..., startFrame] = m_motion1[0, ..., startFrame]
    for (unsigned int i = 0; i < startFrame + 2; i++)
    {
        m_motion3.m_keyFrames[i]->Clone(*(m_motion1.m_keyFrames[i]));
    }
	
    // Estimate root translation and rotation for m_motion3[startFrame + interpFrame + 1]
    // Estimate the speed of the roots during the start and end frames.
    // Find the position and orientation of the end frame roots repectively.
    // TODO: Start your code here...   

	vec3 stVelocity;
	vec3 edVelocity; 
	
	stVelocity = m_motion1.m_keyFrames[startFrame]->GetRootTranslation() - m_motion1.m_keyFrames[startFrame-1]->GetRootTranslation();
	edVelocity = -m_motion2.m_keyFrames[endFrame - 1]->GetRootTranslation() + m_motion2.m_keyFrames[endFrame]->GetRootTranslation();
	vec3 translate = m_motion1.m_keyFrames[startFrame]->GetRootTranslation();

	//make conditional judgement 
	if(((stVelocity[0]>0)&&(edVelocity[0]<0))||((stVelocity[0]<0)&&(edVelocity[0]>0)))
	{
		edVelocity[0] = -1*stVelocity[0]; 
	}
	if(((stVelocity[2]>0)&&(edVelocity[2]<0))||((stVelocity[2]<0)&&(edVelocity[2]>0)))
	{
		edVelocity[2] = -1*stVelocity[2];  
	}
	for(float i = 0; i <= 1; i += (float)1/(interpFrame+1))
	{
		translate = translate +  stVelocity + i * (edVelocity - stVelocity);																	
	}

	translate[1] = m_motion2.GetFrameAt(endFrame)->GetRootTranslation()[1]; 

	Transform newTransF;
	Transform oldTransF;
	Transform convert;

	oldTransF.m_rotation = m_motion2.GetFrameAt(endFrame)->GetJointRotation(0);
	oldTransF.m_translation = m_motion2.GetFrameAt(endFrame)->GetRootTranslation();
	newTransF.m_rotation = m_motion1.GetFrameAt(startFrame)->GetJointRotation(0);
	newTransF.m_translation = translate;

	convert = newTransF * oldTransF.Inverse();																						// Computing the factor to estimate translation and rotation
																															// for motion3
	Transform t;
	Transform tOld;
	Transform tNew;
	t.m_rotation = convert.m_rotation;
	t.m_translation = convert.m_translation;

    // m_motion3[startFrame + interpFrame + 1, ..., totalFrame] = m_motion2[endFrame, ..., frameCount] with proper root translation and rotation
    // TODO: Start your code here...  

	for(unsigned int i = startFrame + interpFrame + 1; i < totalFrame; i++)
	{

		 m_motion3.m_keyFrames[i]->Clone(*(m_motion2.m_keyFrames[endFrame + i-startFrame - interpFrame -1]));
		
		tOld.m_translation = m_motion2.m_keyFrames[endFrame + i-startFrame - interpFrame -1]->GetRootTranslation();
		tOld.m_rotation = m_motion2.m_keyFrames[endFrame + i-startFrame - interpFrame -1]->GetJointRotation(0);

		tNew = t * tOld;																						
		m_motion3.m_keyFrames[i]->SetRootTranslation(vec3(tNew.m_translation[0], tOld.m_translation[1], tNew.m_translation[2])); // Restricting translation in y axis
		m_motion3.m_keyFrames[i]->SetJointRotation(0, tNew.m_rotation);
	} 

    // m_motion3[startFrame + 1, ..., startFrame + interpFrame] are computed from proper interpolations including
    // root translation interpolation
    // joint rotation interpolation using Cubic, SLERP and SQUAD
    float fPerc;
    Frame s0, s1;
    s0.SetJointCount(m_motion3.m_jointCount);
    s1.SetJointCount(m_motion3.m_jointCount);
	Frame::IntermediateFrame(*m_motion3.GetFrameAt(startFrame - 1), *m_motion3.GetFrameAt(startFrame), *m_motion3.GetFrameAt(startFrame + interpFrame + 1), s0);
	Frame::IntermediateFrame(*m_motion3.GetFrameAt(startFrame), *m_motion3.GetFrameAt(startFrame + interpFrame + 1), *m_motion3.GetFrameAt(startFrame + interpFrame + 2), s1);
			
    Frame *pFrame;
	int j = 0; 
    for (int i = startFrame + 1, j = 0; i < startFrame + interpFrame + 1; i++, j++)
    {
        fPerc = float(j + 1) / float(interpFrame + 1);
        pFrame = m_motion3.GetFrameAt(i);
        switch(type)
		{
        case eCUBIC:
            // check the Frame class and the methods inside might be helpful.
            // TODO: cubic interpolation
			pFrame->Cubic(*m_motion3.GetFrameAt(startFrame - 1), *m_motion3.GetFrameAt(startFrame), *m_motion3.GetFrameAt(startFrame + interpFrame + 1), *m_motion3.GetFrameAt(startFrame + interpFrame + 2), *pFrame, fPerc);
            break;
        case eSLERP:
            // TODO: SLERP
			pFrame->Slerp(*m_motion3.GetFrameAt(startFrame), *m_motion3.GetFrameAt(startFrame + interpFrame + 2), *pFrame, fPerc);
            break;
        case eSQUAD:
			pFrame->Squad(*m_motion3.GetFrameAt(startFrame), s0, s1, *m_motion3.GetFrameAt(startFrame + interpFrame + 2), *pFrame, fPerc);
			// TODO: SQUAD
            break;
        }
    }
}





