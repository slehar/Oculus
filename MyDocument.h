//
//  MyDocument.h
//  Oculus
//
//  Created on 9/18/04.
//  Copyright __MyCompanyName__ 2004 . All rights reserved.
//

#import <Cocoa/Cocoa.h>

typedef struct {
	float min;
	float max;
} dynRange;

@interface MyDocument : NSDocument
{
    IBOutlet id sourceView, outView, PSFView;
    IBOutlet id RGBCheckbox;
    IBOutlet id ints, floats1, floats2;
	IBOutlet id processButton;
	IBOutlet id scaleField;
    IBOutlet id scaleSlider;
	IBOutlet id dynRadios1, dynRadios2;
	IBOutlet id dynSlider;
	IBOutlet id dynMinField, dynMaxField;
	IBOutlet id compressionView, compressionViewSourceView, compressionViewOutView;
	IBOutlet id compressionViewSlider, compressionViewField;
	IBOutlet id int3PopUp;
	IBOutlet id sourceSizeField;
	
	NSImage *sourceImage, *outImage, *PSFProcessImage, *PSFImage;
	NSSize imageSize;
	dynRange *dynRanges;
	NSArray *controlNames;
}

- (IBAction)setScale:(id)sender;
- (IBAction)doProcessing:(id)sender;
- (IBAction)setDyn:(id)sender;
- (IBAction)setDynMin:(id)sender;
- (IBAction)setDynMax:(id)sender;
- (IBAction)changeDynRadio:(id)sender;
- (IBAction)changeFloat:(id)sender;
- (IBAction)setControlNames:(id)sender;
- (IBAction)setCompression:(id)sender;

- (void)updateDynSlider;
- (void)scaleBy:(float)percent;
- (void)scalePSFImage;
- (void)setSourceImage:(NSImage *)img;
- (void)setOutImage:(NSImage *)img;
- (void)setPSFImage:(NSImage *)img;
- (void)setPSFProcessImage:(NSImage *)img;
- (void)initDynRanges;

- (NSImage *)outImage;
- (NSView *)outView;
- (NSView *)sourceView;

@end
