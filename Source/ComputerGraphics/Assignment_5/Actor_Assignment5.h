// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Hw5_Shape.h"
#include "Camera/CameraComponent.h"
#include "Components/BillboardComponent.h"
#include "Actor_Assignment5.generated.h"

//���ص���Ϣ
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

	void InitialScene();//��ʼ����Ļ��������ڷŵ�λ

	//���н��й���ɨ�裬���� Texture
	void Render();
	// ������պͷ�������ɫ
	FVector castRay(const FVector& orig, const FVector& dir, int depth);
	// �жϹ����������Ƿ��ཻ
	TOptional<FHit_payload> trace(const FVector& orig, const FVector& dir);

	FVector reflect(const FVector& I, const FVector& N); //����
	FVector refract(const FVector& I, const FVector& N, const float& ior); //����
	float fresnel(const FVector& I, const FVector& N, const float& ior);//���������ϵ��

	void CreateTexture();//���� Texture

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
	float epsilon = 0.001; //����ƫ����Ҫ�Ŵ�ʮ�������������

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
