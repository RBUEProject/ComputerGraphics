// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Hw5_Shape.h"
#include "Camera/CameraComponent.h"
#include "Components/BillboardComponent.h"
#include "Actor_Assignment5.generated.h"

//负载的信息
struct FHit_payload
{
	float tNear;
	uint32 index;
	FVector2D uv;
	AHw5_Shape* hit_obj;
};
UCLASS()
class COMPUTERGRAPHICS_API AActor_Assignment5 : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AActor_Assignment5();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	void InitialScene();//初始化屏幕，把物体摆放到位

	//按行进行光线扫描，更新 Texture
	void Render();
	// 计算光照和反射后的颜色
	FVector castRay(const FVector& orig, const FVector& dir, int depth);
	// 判断光线与物体是否相交
	TOptional<FHit_payload> trace(const FVector& orig, const FVector& dir);

	FVector reflect(const FVector& I, const FVector& N); //反射
	FVector refract(const FVector& I, const FVector& N, const float& ior); //折射
	float fresnel(const FVector& I, const FVector& N, const float& ior);//计算菲涅尔系数

	void CreateTexture();//创建 Texture

	void TextureFromImgArr(const TArray<FColor>& SrcData);

public:
	UPROPERTY(EditAnywhere)
		int32 width = 128;
	UPROPERTY(EditAnywhere)
		int32 height = 96;

	float fov = 90;
	FVector backgroundColor = FVector(0.235294, 0.67451, 0.843137);

	UPROPERTY(EditAnywhere)
		int maxDepth = 5;
	float epsilon = 0.001; //本处偏移需要放大十倍，否则有噪点

	TArray<FVector> lights_pos;
	float lights_Intensity;

	UPROPERTY()
		TArray<AHw5_Shape*> objects;

	UPROPERTY(VisibleAnywhere)
		USceneComponent* rootComp;
	UPROPERTY(VisibleAnywhere)
		UCameraComponent* cameraComp;
	UPROPERTY(VisibleAnywhere)
		UBillboardComponent* light1;
	UPROPERTY(VisibleAnywhere)
		UBillboardComponent* light2;
	UPROPERTY(VisibleAnywhere, BlueprintReadWrite)
		UTexture2D* T_Result;

	bool bFinishedScan = false;
	int32 rowNumber = 0;

	UPROPERTY()
		TArray<FColor> framebuffer;

};
