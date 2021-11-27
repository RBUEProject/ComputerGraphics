// Fill out your copyright notice in the Description page of Project Settings.


#include "Actor_assignment1.h"
#include "Kismet/KismetSystemLibrary.h"
#include "Kismet/KismetMathLibrary.h"
#include "GameFramework/HUD.h"
#include "Kismet/GameplayStatics.h"
#include "MyHUD.h"

//#include "MyHUD.h"
// Sets default values
AActor_assignment1::AActor_assignment1()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	root = CreateDefaultSubobject<USceneComponent>("root");
	SetRootComponent(root);


	//定义的三角形的三个点
	Points.Add(FVector(2.0f,0,-2.0f));
	Points.Add(FVector(0,2.0f,-2.0f));
	Points.Add(FVector(-2.0f,0,-2.0f));
	Triangles.Add(2);
	Triangles.Add(1);
	Triangles.Add(0);

	//初始化
	width = height = 0;
	eye_loc = FVector(0,0,5);
	angle = 0;

	modelMatrix.SetIdentity();
	viewMatrix.SetIdentity();
	projectionMatrix.SetIdentity();

	ponitA.SetTranslation(Points[0]);
	ponitB.SetTranslation(Points[1]);
	ponitC.SetTranslation(Points[2]);
}

// Called when the game starts or when spawned
void AActor_assignment1::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AActor_assignment1::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	FVector2D ViewportSize;
	GetWorld()->GetGameViewport()->GetViewportSize(ViewportSize);
	width = height = ViewportSize.X / 2;

	//绕任意轴旋转MVP变换
	modelMatrix = get_model_matrix_anyAxis(FVector(0, 0, 5), angle);
	viewMatrix = get_view_matrix(eye_loc);
	projectionMatrix = get_projection_matrix(45, 1, 0.1, 50);

	RasterizerDraw();//绘制
	angle+=0.5;
	root->SetWorldRotation(FRotator(0,angle,0));//物体旋转
}

//绘制三角形
void AActor_assignment1::DrawTriangleIn3D()
{
	UKismetSystemLibrary::DrawDebugLine(GetWorld(), Points[0], Points[1], FLinearColor::Green, 0.02f, .2f);
	UKismetSystemLibrary::DrawDebugLine(GetWorld(), Points[1], Points[2], FLinearColor::Green, 0.02f, .2f);
	UKismetSystemLibrary::DrawDebugLine(GetWorld(), Points[2], Points[0], FLinearColor::Green, 0.02f, .2f);
}

FMatrix AActor_assignment1::get_model_matrix_anyAxis(FVector axis, float rotation_angle)
{
	FMatrix model = FMatrix::Identity;//用单位矩阵初始化
	axis.Normalize(0.0001);
	FMatrix N = FMatrix(
		FPlane(0, -axis.Z, axis.Y, 0),
		FPlane(axis.Z, 0, -axis.X, 0),
		FPlane(-axis.Y, axis.X, 0, 0),
		FPlane(0, 0, 0, 0));

	FMatrix rotate4f = FMatrix::Identity * UKismetMathLibrary::DegCos(rotation_angle);
	// nnt = axis x axis的转置
	FMatrix nnT = FMatrix(
		FPlane(axis.X*axis.X, axis.X*axis.Y, axis.X*axis.Z, 0),
		FPlane(axis.Y*axis.X, axis.Y*axis.Y, axis.Y*axis.Z, 0),
		FPlane(axis.Z*axis.X, axis.Z*axis.Y, axis.Z*axis.Z, 0),
		FPlane(0, 0, 0, 0));

	rotate4f += nnT * (1 - UKismetMathLibrary::DegCos(rotation_angle));//DegCos把角度转换成弧度后cos

	rotate4f += N * UKismetMathLibrary::DegSin(rotation_angle);

	rotate4f.M[3][3] = 1;
	model = rotate4f * model;
	return model;
}

FMatrix AActor_assignment1::get_view_matrix(FVector eye_pos)
{
	FMatrix view = FMatrix::Identity;

	FMatrix translate = FMatrix(
		FPlane(1, 0, 0, -eye_pos.X),
		FPlane(0, 1, 0, -eye_pos.Y),
		FPlane(0, 0, 1, -eye_pos.Z),
		FPlane(0, 0, 0, 1));

	view = translate * view;
	return view;
}


FMatrix AActor_assignment1::get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
	FMatrix projection = FMatrix::Identity;

	float t = zNear * UKismetMathLibrary::DegTan(eye_fov / 2);
	float b = -t;
	float r = t * aspect_ratio;
	float l = -r;

	FMatrix translate = FMatrix(
		FPlane(2 * zNear / (r - l), 0, -(r + l) / (r - l), 0),
		FPlane(0, 2 * zNear / (t - b), -(t + b) / (t - b), 0),
		FPlane(0, 0, -(zNear + zFar) / (zNear - zFar), 2 * zNear * zFar / (zNear - zFar)),
		FPlane(0, 0, 1, 0));
	projection = translate * projection;

	return projection;
}

void AActor_assignment1::RasterizerDraw()
{
	FMatrix mvp = projectionMatrix * viewMatrix * modelMatrix;//MVP矩阵
	//f1 f2，用来分割像素点
	float f1 = (100 - 0.1) / 2.0;
	float f2 = (100 + 0.1) / 2.0;
	TArray<FVector4>v;
	//对每个点MVP变换
	for (FVector& p : Points) {
		v.Add(mvp.GetTransposed().TransformFVector4(FVector4(p.X, p.Y, p.Z, 1.0f)));
	}
	//顶点位置
	for (FVector4& vert : v) {
		vert *= 1 / vert.W;
		vert.X = 0.5 * width * (vert.X + 1.0);
		vert.Y = 0.5 * height * (vert.Y + 1.0);
		vert.Z = vert.Z * f1 + f2;
	}
	TArray<FVector> triangleVerts;
	for (FVector4& vert : v) {
		triangleVerts.Add(UKismetMathLibrary::Conv_Vector4ToVector(vert));
	}

	// 调用AHUD 屏幕绘制函数
	AMyHUD* myHUD = Cast<AMyHUD>(UGameplayStatics::GetPlayerController(GetWorld(), 0)->GetHUD());
	if (myHUD) {
		myHUD->rasterize_wireframe(triangleVerts);
	}
}

