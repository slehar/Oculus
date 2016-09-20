//
//  ArrayEditor.m
//  Oculus
//
//  Created by Stephen Poprocki on 6/5/05.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

#import "ArrayEditor.h"


@implementation ArrayEditor

static ArrayEditor *sharedInstance = nil;

+ (ArrayEditor *)sharedInstance
{
    return sharedInstance ? sharedInstance : [[self alloc] init];
}

- (id)init
{
    if(sharedInstance)
    {
        // We just have one instance of the ArrayEditor class, return that one instead
        [self release];
    }
    else if(self = [super init])
    {		
        sharedInstance = self;
    }
    return sharedInstance;
}


- (IBAction)showArrayEditor:(id)sender
{
	[window makeKeyAndOrderFront:self];
}

- (IBAction)setTransparency:(id)sender
{
	[window setAlphaValue:[sender floatValue]];
}

- (char *)arrayWithWidth:(int)w height:(int)h
{
	NSImage *newImg = [[NSImage alloc] initWithSize:NSMakeSize(w,h)];
	[newImg lockFocus];
	[[NSColor grayColor] set];
	NSRect rect = {0, 0, [newImg size].width, [newImg size].height};
	NSRectFill(rect);
	[newImg unlockFocus];
	
	[self setImage:newImg];
	[newImg release];
	bitmap = [NSBitmapImageRep imageRepWithData:[[self image] TIFFRepresentation]];
	return (char *)[bitmap bitmapData];
}

- (void)updateImage
{
	[self setImage:[[[NSImage alloc] initWithData:[bitmap TIFFRepresentation]] autorelease]];
}

- (NSImage *)image { return image; }

- (void)setImage:(NSImage *)newImg
{
	[newImg retain];
	[image release];
	image = newImg;
	
	[arrayView setImage:image];
	[arrayView setFrameSize:NSMakeSize([image size].width, [image size].height)];
	[arrayView setNeedsDisplay:YES];
}

- (NSBitmapImageRep *)bitmap { return bitmap; }

@end
