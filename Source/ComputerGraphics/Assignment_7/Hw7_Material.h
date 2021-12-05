// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Hw7_Global.h"

// ��������
enum  EMaterialType { DIFFUSE };

// ����
class Material {
public:
	EMaterialType m_type;
	//FVector m_color;
	FVector m_emission;
	float ior;
	FVector Kd, Ks;
	float specularExponent;
	//Texture tex;

	Material(EMaterialType t=EMaterialType::DIFFUSE, FVector e=FVector::ZeroVector);
	

public:

	inline EMaterialType getType() { return m_type; }
	//inline FVector getColor() { return m_color; }
	inline FVector getColorAt(double u, double v) { return FVector(); }
	FVector getEmission() {return m_emission;}
	bool hasEmission() {
		if (m_emission.Size() > EPSILON) 
			return true;
		else 
			return false;
	}

	 // sample a ray by Material properties
	// ���ոò��ʵ����ʣ��������䷽���뷨��������ĳ�ֲַ�����һ�����䷽��
	FVector sample(const FVector& wi, const FVector& N);

	// given a ray, calculate the PdF of this ray 
	// ����һ�����䡢���䷽���뷨����������sample �����õ��ó��䷽��ĸ����ܶ�
	float pdf(const FVector& wi, const FVector& wo, const FVector& N);
	
	// given a ray, calculate the contribution of this ray
	// ����һ�����䡢���䷽���뷨������������������µ�f_r ֵ
	FVector eval(const FVector& wi, const FVector& wo, const FVector& N);

	FVector reflect(const FVector& I, const FVector& N); //����
	FVector refract(const FVector& I, const FVector& N); //����
	float fresnel(const FVector& I, const FVector& N); //������
	FVector toWorld(const FVector& a, const FVector& N);
};