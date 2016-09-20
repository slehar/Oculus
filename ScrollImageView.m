//
//  ScrollImageView.m
//  Oculus
//
//  Created on 9/24/04.
//  Copyright 2004 __MyCompanyName__. All rights reserved.
//

#import "ScrollImageView.h"
#import "MyDocument.h"

@implementation ScrollImageView

- (id)initWithFrame:(NSRect)frame
{
    self = [super initWithFrame:frame];
    if(self)
	{
		[self setImageScaling:NSScaleToFit];
    }
    return self;
}

- (void)awakeFromNib
{
	handImage = [NSImage imageNamed:@"hand"];
	handCursor = [[NSCursor alloc] initWithImage:handImage hotSpot:NSMakePoint(7, 7)];
	handGrabImage = [NSImage imageNamed:@"grab"];
	handGrabCursor = [[NSCursor alloc] initWithImage:handGrabImage hotSpot:NSMakePoint(8, 9)];
	[[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(boundsDidChange:) name:NSViewBoundsDidChangeNotification object:[[self enclosingScrollView] contentView]];
}

- (void)boundsDidChange:(NSEvent *)theEvent
{
	NSPoint pt = [[[self enclosingScrollView] contentView] documentVisibleRect].origin;
	[[controller outView] scrollPoint:pt];
	[[controller sourceView] scrollPoint:pt];
}

- (void)resetCursorRects
{
	if(downPt.x==0 && downPt.y==0)
		[self addCursorRect:[self visibleRect] cursor:handCursor];
	else
		[self addCursorRect:[self visibleRect] cursor:handGrabCursor];
}

- (BOOL)acceptsFirstResponder
{
	return YES;
}

- (BOOL)becomeFirstResponder
{
	return YES;
}

- (BOOL)resignFirstResponder
{
	return YES;
}

- (void)mouseDown:(NSEvent *)theEvent
{
	downPt = [theEvent locationInWindow];
	[handGrabCursor set];
	[[self window] invalidateCursorRectsForView:self];
}

- (void)mouseUp:(NSEvent *)theEvent
{
	[handCursor set];
	downPt.x = downPt.y = 0;
	[[self window] invalidateCursorRectsForView:self];
}

- (void)mouseDragged:(NSEvent *)theEvent
{
	NSPoint pt = [theEvent locationInWindow];
	NSPoint dPt;
	NSPoint curPt = [[[self enclosingScrollView] contentView] documentVisibleRect].origin;
	dPt.x = pt.x - downPt.x;
	dPt.y = pt.y - downPt.y;
	[self scrollPoint:NSMakePoint(curPt.x - dPt.x, curPt.y - dPt.y)];
	downPt = pt;
}

@end
