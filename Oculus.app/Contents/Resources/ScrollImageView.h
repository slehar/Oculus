//
//  ScrollImageView.h
//  Oculus
//
//  Created on 9/24/04.
//  Copyright 2004 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface ScrollImageView : NSImageView
{
	IBOutlet id controller;
	NSPoint downPt;
	NSCursor *handCursor;
    NSImage *handImage;
    NSCursor *handGrabCursor;
    NSImage *handGrabImage;
}

@end
