//
//  ArrayEditor.h
//  Oculus
//
//  Created by Stephen Poprocki on 6/5/05.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface ArrayEditor : NSObject {
	IBOutlet id window;
	IBOutlet id arrayView;
	NSImage *image;
	NSBitmapImageRep *bitmap;
}

+ (ArrayEditor *)sharedInstance;

- (IBAction)showArrayEditor:(id)sender;
- (IBAction)setTransparency:(id)sender;
- (char *)arrayWithWidth:(int)w height:(int)h;
- (void)updateImage;
- (NSImage *)image;
- (void)setImage:(NSImage *)newImg;
- (NSBitmapImageRep *)bitmap;

@end
