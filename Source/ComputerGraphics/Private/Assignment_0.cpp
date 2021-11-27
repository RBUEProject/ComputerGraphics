// Fill out your copyright notice in the Description page of Project Settings.


#include "Assignment_0.h"
#include "Kismet/KismetMathLibrary.h"
#include "Kismet/KismetStringLibrary.h"
// Sets default values
AAssignment_0::AAssignment_0()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

}

// Called when the game starts or when spawned
void AAssignment_0::BeginPlay()
{
	Super::BeginPlay();
	float fcos = UKismetMathLibrary::DegCos(45);
	float fsin = UKismetMathLibrary::DegSin(45);
	FPlane row1 = FPlane(fcos,-fsin,1,0);//四元数处理二维矩阵变换
	FPlane row2 = FPlane(fsin,fcos,2,0);
	FPlane row3 = FPlane(0,0,1,0);
	FPlane row4 = FPlane(0,0,0,0);
	FMatrix TransMatrix = FMatrix(row1,row2,row3,row4);

	FVector4 originPos = FVector4(2,1,1,0);
	TransMatrix = TransMatrix.GetTransposed();//求矩阵的转置
	FVector4 res = TransMatrix.TransformFVector4(originPos);
	//TransMatrix(T)*originPos == TransMatrix * originPos(T) 整体再转置

	UE_LOG(LogTemp, Warning, TEXT("%s"), *TransMatrix.ToString());
	UE_LOG(LogTemp, Warning, TEXT("[ %.2f, %.2f, %.2f ]"), res.X, res.Y, res.Z);
}

// Called every frame
void AAssignment_0::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
}

