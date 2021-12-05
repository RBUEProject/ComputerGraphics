// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Hw7_Global.h"
#include "Hw7_Object.h"
#include "BVHAccelTree.h"
#include "Hw7_Material.h"
#include "Hw7_Triangle.generated.h"

UCLASS()
class COMPUTERGRAPHICS_API UHw7_Triangle : public UObject, public IHw7_Object
{
	GENERATED_BODY()
public:
	FVector v0, v1, v2; // vertices A, B ,C , counter-clockwise order
	FVector e1, e2;     // 2 edges v1-v0, v2-v0;
	FVector t0, t1, t2; // texture coords
	FVector normal;
	float area;
	Material* m;
	Bounds3 bounding_box;

	UHw7_Triangle() {}

	void Init(FVector _v0, FVector _v1, FVector _v2, Material* _m = nullptr) {
		v0 = _v0;
		v1 = _v1;
		v2 = _v2;
		e1 = v1 - v0;
		e2 = v2 - v0;
		normal = FVector::CrossProduct(e2, e1); //ue4 ���������⣬��Ҫ����λ��
		area = normal.Size() * 0.5f;
		normal.Normalize();
		m = _m;
		bounding_box = Bounds3(v0, v1);
		bounding_box.Union(v2);
	}

	bool intersect(const Ray& ray) override { return true; }
	bool intersect(const Ray& ray, float& tnear, uint32_t& index) const override { return true; }

	Intersection getIntersection(Ray ray) override;

	void getSurfaceProperties(
		const FVector& P, const FVector& I,
		const uint32_t& index, const FVector2D& uv,
		FVector& N, FVector2D& st) const override {
		N = normal;
	};

	FVector evalDiffuseColor(const FVector2D&) const override { return FVector(0.5, 0.5, 0.5); }

	Bounds3 getBounds() const override { return bounding_box; }
	
	//�����̡߳�����
	void Sample(Intersection& pos, float& pdf) {
		float x = std::sqrt(UKismetMathLibrary::RandomFloat()), y = UKismetMathLibrary::RandomFloat();
		pos.coords = v0 * (1.0f - x) + v1 * (x * (1.0f - y)) + v2 * (x * y);
		pos.normal = this->normal;
		pdf = 1.0f / area;
	}
	float getArea() {return area;}
	bool hasEmit() {return m->hasEmission();}
};


UCLASS()
class COMPUTERGRAPHICS_API AHw7_MeshTriangle : public AActor, public IHw7_Object
{
	GENERATED_BODY()

public:
	AHw7_MeshTriangle();

public:

	Bounds3 bounding_box;

	UPROPERTY()
		TArray<FVector> vertices;
	UPROPERTY()
		TArray<uint32> vertexIndex;
	UPROPERTY()
		TArray<FVector2D> stCoordinates;

	UPROPERTY()
		TArray<TScriptInterface<IHw7_Object>> triangles;

	UPROPERTY()
		UBVHAccelTree* bvhTree;
	float area = 0;
	Material* m;

	UPROPERTY(VisibleAnywhere)
		UStaticMeshComponent* meshComp;

public:
	void Init(Material* _m = nullptr);

	void getStaticMeshData();

	//�� mesh ��ֳ������湹�� BVH Tree
	void buildBVHTree();

	// �жϹ����Ƿ���������ཻ
	bool rayTriangleIntersect(
		const FVector& v0, const FVector& v1,
		const FVector& v2, const FVector& orig,
		const FVector& dir, float& tnear,
		float& u, float& v) const;


	bool intersect(const Ray& ray) override { return true; }

	// �ж��Ƿ��ཻ
	bool intersect(const Ray& ray, float& tnear, uint32_t& index) const override;

	// ��ȡ�ռ�����
	Bounds3 getBounds() const override { return bounding_box; }

	// ��ȡ������Ϣ
	void getSurfaceProperties(
		const FVector& P, const FVector& I,
		const uint32_t& index, const FVector2D& uv,
		FVector& N, FVector2D& st) const override;

	// ��ȡ���Ӧ����ɫ
	FVector evalDiffuseColor(const FVector2D& st) const override;

	// ��ȡ�ཻ��Ϣ
	Intersection getIntersection(Ray ray) override;

	//�����̡߳�����
	void Sample(Intersection& pos, float& pdf) {
		bvhTree->Sample(pos, pdf);
		pos.emit = m->getEmission();
	}
	float getArea() { return area;}
	bool hasEmit() {return m->hasEmission();}



};
