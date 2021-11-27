// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Actor_assignment1.generated.h"

UCLASS()
class COMPUTERGRAPHICS_API AActor_assignment1 : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AActor_assignment1();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;
	void DrawTriangleIn3D();

	FMatrix get_view_matrix(FVector eye_pos);//视口变换
	FMatrix get_model_matrix_anyAxis(FVector axis, float rotation_angle);//绕任意轴模型变换
	FMatrix get_projection_matrix(float eye_pos,float aspect_ratio,float zNear,float zFar);//投影变换

	void RasterizerDraw();

public:
		FTransform ponitA;
		FTransform ponitB;
		FTransform ponitC;

	UPROPERTY()
		USceneComponent* root;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (MakeEditWidget))
		TArray<FVector> Points;
	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (MakeEditWidget))
		TArray<int32> Triangles;

	int32 width, height;
	FVector eye_loc;
	float angle;

	FMatrix modelMatrix;
	FMatrix viewMatrix;
	FMatrix projectionMatrix;
};
