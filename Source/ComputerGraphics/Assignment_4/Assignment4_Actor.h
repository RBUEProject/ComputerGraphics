// Fill out your copyright notice in the Description page of Project Settings.

#pragma once
#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/SceneComponent.h"
#include "Assignment4_Actor.generated.h"

UCLASS()
class COMPUTERGRAPHICS_API AAssignment4_Actor : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AAssignment4_Actor();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;
	
	UFUNCTION(BlueprintCallable)
	void naive_bezier();
	UFUNCTION(BlueprintCallable)
	void bezier();
	UFUNCTION(BlueprintCallable)
	FVector recursive_bezier(TArray<FVector>&points,float t);//µÝ¹ébezier

public:
	UPROPERTY(VisibleAnywhere)
	USceneComponent* root;
	UPROPERTY(VisibleAnywhere, Category = "Points")
		UStaticMeshComponent* point0;
	UPROPERTY(VisibleAnywhere, Category = "Points")
		UStaticMeshComponent* point1;
	UPROPERTY(VisibleAnywhere, Category = "Points")
		UStaticMeshComponent* point2;
	UPROPERTY(VisibleAnywhere, Category = "Points")
		UStaticMeshComponent* point3;
	UPROPERTY(VisibleAnywhere, Category = "Points")
		UStaticMeshComponent* point4;

	UPROPERTY();
	TArray<FVector> m_points;

	UPROPERTY(EditAnywhere);
	bool m_bUseRecursiveBezier;

};
