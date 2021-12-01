// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "Slate/SlateBrushAsset.h"
#include "Components/Image.h"
#include "Components/Button.h"
#include "Rasterizer2Widget.generated.h"

/**
 * 
 */
UCLASS()
class COMPUTERGRAPHICS_API URasterizer2Widget : public UUserWidget
{
	GENERATED_BODY()
protected:
	virtual void NativePreConstruct();
public:
	//UFUNCTION(BlueprintCallable, Category = "GAMES101 | 2")
	//	void DrawBackground(FPaintContext& Context) const;

	//void Draw(TArray<FVector>& frame_buf);

	UFUNCTION()
		void SwitchMASS();


public:
	UPROPERTY()
		USlateBrushAsset* SlateBrushAsset;

	UPROPERTY(EditAnywhere, meta = (BindWidget))
		UImage* TCanvas;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (BindWidget))
		UButton* Btn_Mass;

private:
	FVector2D CanvasSize;
	bool bDrawSwitch = false;
	int32 drawInterval = 10;

};
