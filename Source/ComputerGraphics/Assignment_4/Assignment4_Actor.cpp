// Fill out your copyright notice in the Description page of Project Settings.


#include "Assignment4_Actor.h"
#include "DrawDebugHelpers.h"
#include "cmath"
#include "algorithm"
#include "Kismet/KismetSystemLibrary.h"
// Sets default values
AAssignment4_Actor::AAssignment4_Actor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	root = CreateDefaultSubobject<USceneComponent>("root");
	SetRootComponent(root);

	point0 = CreateDefaultSubobject<UStaticMeshComponent>("point0");
	point1 = CreateDefaultSubobject<UStaticMeshComponent>("point1");
	point2 = CreateDefaultSubobject<UStaticMeshComponent>("point2");
	point3 = CreateDefaultSubobject<UStaticMeshComponent>("point3");
	point4 = CreateDefaultSubobject<UStaticMeshComponent>("point4");

	point0->SetupAttachment(root);
	point1->SetupAttachment(root);
	point2->SetupAttachment(root);
	point3->SetupAttachment(root);
	point4->SetupAttachment(root);

	m_points.Init(FVector::ZeroVector,5);//用5个0向量初始化
	m_bUseRecursiveBezier = false;
}

// Called when the game starts or when spawned
void AAssignment4_Actor::BeginPlay()
{
	Super::BeginPlay();
	m_points[0] = point0->GetComponentLocation();
	m_points[1] = point1->GetComponentLocation();
	m_points[2] = point2->GetComponentLocation();
	m_points[3] = point3->GetComponentLocation();
	m_points[4] = point4->GetComponentLocation();

	if (!m_bUseRecursiveBezier)
		naive_bezier();
	else
		bezier();
}

// Called every frame
void AAssignment4_Actor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void AAssignment4_Actor::naive_bezier()//简单贝塞尔
{
	FVector& p_0 = m_points[0];
	FVector& p_1 = m_points[1];
	FVector& p_2 = m_points[2];
	FVector& p_3 = m_points[3];
	FVector& p_4 = m_points[4];
	for (double t = 0.0; t < 1.0; t += 0.001) {
		auto point = std::pow(1-t,4)*p_0+4*t*std::pow(1-t,3)*p_1
		+6 * std::pow(t, 2) * std::pow((1 - t), 2) * p_2 
		+ 4 * std::pow(t, 3) * (1 - t) * p_3 + std::pow(t, 4) * p_4;
		DrawDebugPoint(GetWorld(), point, 2.0f, FColor::Red, true, 5.0f);
	}
}

void AAssignment4_Actor::bezier()//递归贝塞尔  德卡斯特里奥
{
	for (double t = 0.0; t <= 1.0; t += 0.001)
	{
		FVector point = recursive_bezier(m_points, t);
		DrawDebugPoint(GetWorld(), point, 2.0f, FColor::Green, true, 5.0f);
	}
}

FVector AAssignment4_Actor::recursive_bezier(TArray<FVector>&points, float t)
{
	if (points.Num() < 3) {
		return (1 - t) * points[0] + t * points[1];
	}
	TArray<FVector> newPoint;
	for (int i = 0; i < points.Num() - 1; i++) {
		newPoint.Add((1 - t) * points[i] + t * points[i + 1]);
	}
	return recursive_bezier(newPoint, t);
}

