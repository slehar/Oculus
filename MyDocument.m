//
//  MyDocument.m
//  Oculus
//
//  Created on 9/18/04.
//  Copyright __MyCompanyName__ 2004 . All rights reserved.
//

#import "deblurring applications.h"
#import "MyDocument.h"
#import "ArrayEditor.h"
#import "math.h"

#define kMaxParams 14
#define kPSFDim 150

@implementation MyDocument

- (void)dealloc
{
	[sourceImage release];
	[outImage release];
	[PSFImage release];
	[PSFProcessImage release];
	[controlNames release];
	free(dynRanges);
	
	[super dealloc];
}

- (IBAction)doProcessing:(id)sender
{

	NSImage *resized;
	NSBitmapImageRep *resizedBitmap;
	char *resizedBitmapPointer;
	NSRect resizeRect = {0,0,kPSFDim,kPSFDim};
	NSSize resizeSize = {kPSFDim,kPSFDim};
	NSBitmapImageRep *bitmap = [NSBitmapImageRep imageRepWithData:[sourceImage TIFFRepresentation]];

	unsigned int pixelsWide = [bitmap pixelsWide];
	unsigned int pixelsHigh = [bitmap pixelsHigh];
	unsigned int samplesPerPixel = [bitmap samplesPerPixel];
	char *bitmapPointer = (char *)[bitmap bitmapData];
	char *arrayPointer = [[ArrayEditor sharedInstance] arrayWithWidth:pixelsWide height:pixelsHigh];
    
	DEBLURINPUT deblurInput = {[[ints cellAtRow:0 column:0] intValue],
							   [[ints cellAtRow:1 column:0] intValue],
							   [int3PopUp indexOfSelectedItem],
							   [RGBCheckbox state],
							   [[floats1 cellAtRow:0 column:0] floatValue],
							   [[floats1 cellAtRow:1 column:0] floatValue],
							   [[floats1 cellAtRow:2 column:0] floatValue],
							   [[floats1 cellAtRow:3 column:0] floatValue],
							   [[floats1 cellAtRow:4 column:0] floatValue],
							   [[floats1 cellAtRow:5 column:0] floatValue],
							   [[floats1 cellAtRow:6 column:0] floatValue],
							   [[floats2 cellAtRow:0 column:0] floatValue],
							   [[floats2 cellAtRow:1 column:0] floatValue],
							   [[floats2 cellAtRow:2 column:0] floatValue],
							   [[floats2 cellAtRow:3 column:0] floatValue],
							   [[floats2 cellAtRow:4 column:0] floatValue],
							   [[floats2 cellAtRow:5 column:0] floatValue],
							   [[floats2 cellAtRow:6 column:0] floatValue]};

	processImage(samplesPerPixel, pixelsWide, pixelsHigh, bitmapPointer, arrayPointer, &deblurInput);
	// is this the fastest way to do it? seems like a lot of overhead, but I can't find better.
	[self setOutImage:[[[NSImage alloc] initWithData:[bitmap TIFFRepresentation]] autorelease]];

	[[ArrayEditor sharedInstance] updateImage];

	// set the PSF image. resize the source to resizeSize and process that image to avoid overhead
	resized = [[NSImage alloc] initWithSize:resizeSize];
	[resized lockFocus];
	[bitmap drawInRect:resizeRect];
	resizedBitmap = [[NSBitmapImageRep alloc] initWithFocusedViewRect:resizeRect];
	[resized unlockFocus];
	resizedBitmapPointer = (char *)[resizedBitmap bitmapData];
	pixelsWide = [resizedBitmap pixelsWide];
	pixelsHigh = [resizedBitmap pixelsHigh];
	samplesPerPixel = [resizedBitmap samplesPerPixel];
	deblurInput.inputint3 = 13;
	deblurInput.inputint4 = 1;
	processImage(samplesPerPixel, pixelsWide, pixelsHigh, resizedBitmapPointer, (char *)nil, &deblurInput);
	[self setPSFImage:[[[NSImage alloc] initWithData:[resizedBitmap TIFFRepresentation]] autorelease]];
	[self setPSFProcessImage:[[[NSImage alloc] initWithData:[resizedBitmap TIFFRepresentation]] autorelease]];
	[resizedBitmap release];
	[resized release];
	[self scalePSFImage];

}

- (IBAction)setScale:(id)sender
{
	if([sender isEqual:scaleSlider])
		[scaleField setFloatValue:[sender floatValue]];
	
	[self scaleBy:[sender floatValue]];
	[self updateChangeCount:NSChangeDone];
}

- (void)scaleBy:(float)percent
{
	float factor = percent/100;
	
	[sourceView setFrameSize:NSMakeSize(factor*imageSize.width, factor*imageSize.height)];
	[sourceView setNeedsDisplay:YES];	
	[outView setFrameSize:NSMakeSize(factor*imageSize.width, factor*imageSize.height)];
	[outView setNeedsDisplay:YES];
	[self scalePSFImage];
}

- (void)scalePSFImage
{
	float factor = [scaleSlider floatValue]/100;
	NSBitmapImageRep *PSFBitmap = [[NSBitmapImageRep imageRepWithData:[PSFProcessImage TIFFRepresentation]] copy];
	NSBitmapImageRep *resizedPSFBitmap;
	NSImage *resizedPSF;
	NSRect resizeRect = {0,0,kPSFDim,kPSFDim};
	NSSize resizeSize = {kPSFDim,kPSFDim};
	
	// scale and center-crop the PSF
	resizedPSF = [[NSImage alloc] initWithSize:resizeSize];
	[resizedPSF lockFocus];
	[PSFBitmap drawInRect:NSMakeRect(-factor*kPSFDim/2+kPSFDim/2,-factor*kPSFDim/2+kPSFDim/2,factor*kPSFDim,factor*kPSFDim)];
	resizedPSFBitmap = [[NSBitmapImageRep alloc] initWithFocusedViewRect:resizeRect];
	[resizedPSF unlockFocus];
	[self setPSFImage:[[[NSImage alloc] initWithData:[resizedPSFBitmap TIFFRepresentation]] autorelease]];
	[resizedPSFBitmap release];
	[resizedPSF release];
}

- (id)dynFloatField
{
	if([dynRadios1 selectedRow] != -1)
		return [floats1 cellAtRow:[dynRadios1 selectedRow] column:0];
	else if([dynRadios2 selectedRow] != -1)
		return [floats2 cellAtRow:[dynRadios2 selectedRow] column:0];
}

// the tag is 1-7 for floats1 and 8-14 for floats2
- (int)dynFloatFieldTag
{
	return ([dynRadios1 selectedRow] != -1) ? // is something in the 1st radios selected?
		[[floats1 cellAtRow:[dynRadios1 selectedRow] column:0] tag] :
		[[floats2 cellAtRow:[dynRadios2 selectedRow] column:0] tag] ;
}

- (IBAction)setDyn:(id)sender
{
	[[self dynFloatField] setFloatValue:[sender floatValue]];
	[self doProcessing:nil];
}

- (IBAction)setDynMin:(id)sender
{
	float curF = [[self dynFloatField] floatValue];
	float newF = [sender floatValue];
	if(newF > curF)
		newF = curF;
	[dynSlider setMinValue:newF];
	dynRanges[[self dynFloatFieldTag]-1].min = newF;
	[sender setFloatValue:newF];
}

- (IBAction)setDynMax:(id)sender
{
	float curF = [[self dynFloatField] floatValue];
	float newF = [sender floatValue];
	if(newF < curF)
		newF = curF;
	[dynSlider setMaxValue:newF];
	dynRanges[[self dynFloatFieldTag]-1].max = newF;
	[sender setFloatValue:newF];
}

- (IBAction)changeDynRadio:(id)sender
{
	// toggle the other rodios off
	if([sender isEqual:dynRadios1])
		[dynRadios2 deselectAllCells];
	else
		[dynRadios1 deselectAllCells];
	
	[self updateDynSlider];
}

- (IBAction)changeFloat:(id)sender
{
	[self updateDynSlider];
	[self doProcessing:nil];
}

- (void)updateDynSlider
{
	float newF = [[self dynFloatField] floatValue];
	int row = [self dynFloatFieldTag]-1;
	
	// make sure the dynSlider's range contains the new value
	if(newF > dynRanges[row].max)
	{
		dynRanges[row].max = newF;
	}
	if(newF < dynRanges[row].min)
	{
		dynRanges[row].min = newF;
	}
	[dynSlider setMaxValue:dynRanges[row].max];
	[dynMaxField setFloatValue:dynRanges[row].max];
	[dynSlider setMinValue:dynRanges[row].min];
	[dynMinField setFloatValue:dynRanges[row].min];
	[dynSlider setFloatValue:newF];
}

- (void)setSourceImage:(NSImage *)img
{
	[img retain];
	[sourceImage release];
	sourceImage = img;
	imageSize = [img size];
	
	[sourceView setImage:sourceImage];
	[self updateChangeCount:NSChangeDone];
}

- (void)setOutImage:(NSImage *)img
{
	[img retain];
	[outImage release];
	outImage = img;
	imageSize = [img size];
	
	[outView setImage:outImage];
	[self updateChangeCount:NSChangeDone];
}

- (void)setPSFImage:(NSImage *)img
{
	[img retain];
	[PSFImage release];
	PSFImage = img;
	
	[PSFView setImage:PSFImage];
}

- (void)setPSFProcessImage:(NSImage *)img
{
	[img retain];
	[PSFProcessImage release];
	PSFProcessImage = img;
}

//	Document Code

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"MyDocument";
}

- (NSImage *)loadImage:(id)sender
{
    NSOpenPanel *oPanel = [NSOpenPanel openPanel];
    NSString *aFile;
    int result;
	
    [oPanel setTitle:@"Choose image to load"];
    [oPanel setAllowsMultipleSelection:NO];
    result = [oPanel runModalForDirectory:nil file:nil types:[NSImage imageFileTypes]];
    if(result == NSOKButton)
    {
		NSImage *img;
        aFile = [[oPanel filenames] objectAtIndex:0];
        img = [[NSImage alloc] initWithContentsOfFile:aFile];
		[self setSourceImage:img];
		[self setOutImage:img];
		
		[sourceSizeField setStringValue:[NSString stringWithFormat:@"Source image size: %d x %d pixels",
										(int)[img size].width, (int)[img size].height]];
		
		[img release];
    }
}

- (void)initDynRanges
{
	int i;
	dynRanges = malloc(kMaxParams*sizeof(dynRange));
	for(i=0; i<kMaxParams; i++)
	{
		dynRanges[i].min = 0.0;
		dynRanges[i].max = 2.0;
	}
}

- (void)setInt3PopUp
{
	int i;
	[int3PopUp removeAllItems];
	for(i=0; i<[controlNames count]; i++)
	{
		NSDictionary *functionDict = [controlNames objectAtIndex:i];
		NSEnumerator *functionDictEnum = [functionDict keyEnumerator];
		NSString *functionName = [NSString stringWithFormat:@"%d %@", i, [functionDictEnum nextObject]];
		[int3PopUp addItemWithTitle:functionName];
	}
}

- (void)resetControlNames
{
	int i;
	for(i=1; i<=[ints numberOfRows]; i++)
	{
		[[ints cellAtIndex:i-1] setTitle:[NSString stringWithFormat:@"Int %d:", i]];
	}
	for(i=1; i<=[floats1 numberOfRows]; i++)
	{
		[[floats1 cellAtIndex:i-1] setTitle:[NSString stringWithFormat:@"Float %d:", i]];
	}
	for(i=1; i<=[floats2 numberOfRows]; i++)
	{
		[[floats2 cellAtIndex:i-1] setTitle:[NSString stringWithFormat:@"Float %d:", i+kMaxParams/2]];
	}
}

- (IBAction)setControlNames:(id)sender
{
	NSDictionary *functionNumDict = [controlNames objectAtIndex:[int3PopUp indexOfSelectedItem]];
	NSDictionary *functionDict = [functionNumDict objectForKey:[[functionNumDict allKeys] objectAtIndex:0]];
	NSDictionary *intsDict = [functionDict objectForKey:@"ints"];
	NSDictionary *floatsDict = [functionDict objectForKey:@"floats"];
	NSEnumerator *enumerator;
	NSString *key;
	
	[self resetControlNames];
	
	enumerator = [intsDict keyEnumerator];
	while(key = [enumerator nextObject])
	{
		int intNum = [key intValue];
		id formCell = [ints cellAtIndex:intNum-1];
		[formCell setTitle:[NSString stringWithFormat:@"%@:", [intsDict objectForKey:key]]];
	}
	
	enumerator = [floatsDict keyEnumerator];
	while(key = [enumerator nextObject])
	{
		int floatNum = [key intValue];
		id formCell;
		if(floatNum <= kMaxParams/2)
			formCell = [floats1 cellAtIndex:floatNum-1];
		else
			formCell = [floats2 cellAtIndex:floatNum-(kMaxParams/2+1)];
		[formCell setTitle:[NSString stringWithFormat:@"%@:", [floatsDict objectForKey:key]]];
	}
    
	[self doProcessing:nil];
}

- (void)windowControllerDidLoadNib:(NSWindowController *) aController
{
	NSString *controlNamesPath = [[NSBundle mainBundle] pathForResource:@"controlNames" ofType:@""];
	controlNames = [[NSArray arrayWithContentsOfFile:controlNamesPath] retain];
	
    [super windowControllerDidLoadNib:aController];
	
	// bring up load image panel
	if(!sourceImage)
		[self loadImage:nil];
	
	[self initDynRanges];
	
	[self setScale:scaleSlider]; // scale the default size in the nib
	
	[self setInt3PopUp];
	[self setControlNames:nil];
}


- (NSData *)dataRepresentationOfType:(NSString *)aType
{
	if([aType isEqualToString:@"JPEG Image"])
	{
		NSDictionary *properties = [NSDictionary dictionaryWithObjectsAndKeys:
			[NSNumber numberWithFloat:[compressionViewSlider floatValue]], NSImageCompressionFactor, nil];
		
		NSBitmapImageRep *imgRep = [NSBitmapImageRep imageRepWithData:[outImage TIFFRepresentation]];
		return [imgRep representationUsingType:NSJPEGFileType properties:properties];
	}
	else if([aType isEqualToString:@"TIFF Document"])
	{
		return [outImage TIFFRepresentation];
	}
	else if([aType isEqualToString:@"Portable Network Graphics Image"])
	{
		NSBitmapImageRep *imgRep = [NSBitmapImageRep imageRepWithData:[outImage TIFFRepresentation]];
		return [imgRep representationUsingType:NSPNGFileType properties:nil];
	}
	/*
	else if([aType isEqualToString:@"Oculus Document"])
	{
		NSMutableDictionary *dictionary = [NSMutableDictionary dictionary];

		[dictionary setObject:[NSNumber numberWithFloat:scale] forKey:@"scale"];
		[dictionary setObject:[NSNumber numberWithInt:RGBProcessing] forKey:@"RGBProcessing"];
		[dictionary setObject:[sourceImage TIFFRepresentation] forKey:@"sourceImage"];
		[dictionary setObject:[outImage TIFFRepresentation] forKey:@"outImage"];

		return [[dictionary description] dataUsingEncoding:NSASCIIStringEncoding];
	}
	*/
	else
		return nil;
}

- (BOOL)prepareSavePanel:(NSSavePanel *)savePanel
{
	NSView *typeView = [[savePanel accessoryView] retain];
	[compressionViewSourceView setFrameSize:[outImage size]];
	[compressionViewSourceView setImage:outImage];
	[compressionViewOutView setFrameSize:[outImage size]];
	[compressionViewOutView setImage:outImage];
	[compressionViewSlider setFloatValue:1.0];
	[compressionViewField setFloatValue:1.0];
	[savePanel setAccessoryView:compressionView];
	[compressionView addSubview:typeView];
	[typeView release];
	return YES;
}

- (IBAction)setCompression:(id)sender
{
	float compressionFactor = [sender floatValue];
	NSImage *compressionViewOutImage;
	NSDictionary *properties = [NSDictionary dictionaryWithObjectsAndKeys:
		[NSNumber numberWithFloat:compressionFactor], NSImageCompressionFactor, nil];
	
	NSBitmapImageRep *imgRep = [NSBitmapImageRep imageRepWithData:[[compressionViewSourceView image] TIFFRepresentation]];
	compressionViewOutImage = [[NSImage alloc] initWithData:[imgRep representationUsingType:NSJPEGFileType properties:properties]];
	[compressionViewOutView setImage:compressionViewOutImage];
	[compressionViewOutImage release];
	
	[compressionViewField setFloatValue:compressionFactor];
}

- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
	/*
	if([aType isEqualToString:@"Oculus Document"])
	{
		NSString *string = [[NSString allocWithZone:[self zone]] initWithData:data encoding:NSASCIIStringEncoding];
		NSDictionary *dictionary = [string propertyList];
		[string release];
		
		scale = [[dictionary objectForKey:@"scale"] floatValue];
		RGBProcessing = [[dictionary objectForKey:@"RGBProcessing"] intValue];
		[self setSourceImage:[[[NSImage alloc] initWithData:[dictionary objectForKey:@"sourceImage"]] autorelease]];
		[self setOutImage:[[[NSImage alloc] initWithData:[dictionary objectForKey:@"outImage"]] autorelease]];
		
		return YES;
	}
	*/
	return NO;
}


- (void)windowWillClose:(NSNotification *)aNotification
{
	[[NSNotificationCenter defaultCenter] removeObserver:outView];
	[[NSNotificationCenter defaultCenter] removeObserver:sourceView];
}

// accessor methods
- (NSImage *)outImage { return outImage; }
- (NSView *)outView { return outView; }
- (NSView *)sourceView { return sourceView; }

@end
