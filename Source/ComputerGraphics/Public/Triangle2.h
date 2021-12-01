// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"

/**
 * 
 */
class COMPUTERGRAPHICS_API Triangle2
{
public:
	Triangle2();
	~Triangle2();

	TArray<FVector> v; 
	FVector color[3]; 
	FVector2D tex_coords[3]; 
	FVector normal[3];

	void setVertex(int ind, FVector ver); 
	void setNormal(int ind, FVector n); 
	void setColor(int ind, float r, float g, float b); 
	FVector getColor() const { return color[0] * 255; }
	void setTexCoord(int ind, float s, float t);
	TArray<FVector4> toVector4() const;
};
