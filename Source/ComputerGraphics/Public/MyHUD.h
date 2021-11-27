// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/HUD.h"
#include "MyHUD.generated.h"

/**
 * 
 */
UCLASS()
class COMPUTERGRAPHICS_API AMyHUD : public AHUD
{
	GENERATED_BODY()
public:
	virtual void DrawHUD() override;
	// �������ڻ���
	void rasterize_wireframe(TArray<FVector>& t);
	// �洢�����ε�������
	TArray<FVector> TriangleVerts;
};
