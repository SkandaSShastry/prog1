/* classes */ 

// Color constructor
class Color {
    constructor(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end try
        
        catch (e) {
            console.log(e);
        }
    } // end Color constructor

        // Color change method
    change(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end throw
        
        catch (e) {
            console.log(e);
        }
    } // end Color change method
} // end color class


/* utility functions */

// draw a pixel at x,y using color
function drawPixel(imagedata,x,y,color) {
    try {
        if ((typeof(x) !== "number") || (typeof(y) !== "number"))
            throw "drawpixel location not a number";
        else if ((x<0) || (y<0) || (x>=imagedata.width) || (y>=imagedata.height))
            throw "drawpixel location outside of image";
        else if (color instanceof Color) {
            var pixelindex = (y*imagedata.width + x) * 4;
            imagedata.data[pixelindex] = color.r;
            imagedata.data[pixelindex+1] = color.g;
            imagedata.data[pixelindex+2] = color.b;
            imagedata.data[pixelindex+3] = color.a;
        } else 
            throw "drawpixel color is not a Color";
    } // end try
    
    catch(e) {
        console.log(e);
    }
} // end drawPixel
    
// draw random pixels
function drawRandPixels(context) {
    var c = new Color(0,0,0,0); // the color at the pixel: black
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.01;
    var numPixels = (w*h)*PIXEL_DENSITY; 
    
    // Loop over 1% of the pixels in the image
    for (var x=0; x<numPixels; x++) {
        c.change(Math.random()*255,Math.random()*255,
            Math.random()*255,255); // rand color
        drawPixel(imagedata,
            Math.floor(Math.random()*w),
            Math.floor(Math.random()*h),
                c);
    } // end for x
    context.putImageData(imagedata, 0, 0);
} // end draw random pixels

// get the input ellipsoids from the standard class URL
function getInputEllipsoids() {
    const INPUT_ELLIPSOIDS_URL = 
        "https://ncsucgclass.github.io/prog1/ellipsoids.json";
        
    // load the ellipsoids file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_ELLIPSOIDS_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input ellipses file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input ellipsoids

//get the input triangles from the standard class URL
function getInputTriangles() {
    const INPUT_TRIANGLES_URL = 
        "https://ncsucgclass.github.io/prog1/triangles2.json";
        
    // load the triangles file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_TRIANGLES_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input triangles file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input triangles

//get the input boxex from the standard class URL
function getInputBoxes() {
    const INPUT_BOXES_URL = 
        "https://ncsucgclass.github.io/prog1/boxes.json";
        
    // load the boxes file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",INPUT_BOXES_URL,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input boxes file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input boxes

// put random points in the ellipsoids from the class github
function drawRandPixelsInInputEllipsoids(context) {
    var inputEllipsoids = getInputEllipsoids();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputEllipsoids != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var cx = 0; var cy = 0; // init center x and y coord
        var ellipsoidXRadius = 0; // init ellipsoid x radius
        var ellipsoidYRadius = 0; // init ellipsoid y radius
        var numEllipsoidPixels = 0; // init num pixels in ellipsoid
        var c = new Color(0,0,0,0); // init the ellipsoid color
        var n = inputEllipsoids.length; // the number of input ellipsoids
        //console.log("number of ellipses: " + n);

        // Loop over the ellipsoids, draw rand pixels in each
        for (var e=0; e<n; e++) {
            cx = w*inputEllipsoids[e].x; // ellipsoid center x
            cy = h*inputEllipsoids[e].y; // ellipsoid center y
            ellipsoidXRadius = Math.round(w*inputEllipsoids[e].a); // x radius
            ellipsoidYRadius = Math.round(h*inputEllipsoids[e].b); // y radius
            numEllipsoidPixels = ellipsoidXRadius*ellipsoidYRadius*Math.PI; // projected ellipsoid area
            numEllipsoidPixels *= PIXEL_DENSITY; // percentage of ellipsoid area to render to pixels
            numEllipsoidPixels = Math.round(numEllipsoidPixels);
            //console.log("ellipsoid x radius: "+ellipsoidXRadius);
            //console.log("ellipsoid y radius: "+ellipsoidYRadius);
            //console.log("num ellipsoid pixels: "+numEllipsoidPixels);
            c.change(
                inputEllipsoids[e].diffuse[0]*255,
                inputEllipsoids[e].diffuse[1]*255,
                inputEllipsoids[e].diffuse[2]*255,
                255); // ellipsoid diffuse color
            for (var p=0; p<numEllipsoidPixels; p++) {
                do {
                    x = Math.random()*2 - 1; // in unit square 
                    y = Math.random()*2 - 1; // in unit square
                } while (Math.sqrt(x*x + y*y) > 1) // a circle is also an ellipse
                drawPixel(imagedata,
                    cx+Math.round(x*ellipsoidXRadius),
                    cy+Math.round(y*ellipsoidYRadius),c);
                //console.log("color: ("+c.r+","+c.g+","+c.b+")");
                //console.log("x: "+Math.round(w*inputEllipsoids[e].x));
                //console.log("y: "+Math.round(h*inputEllipsoids[e].y));
            } // end for pixels in ellipsoid
        } // end for ellipsoids
        context.putImageData(imagedata, 0, 0);
    } // end if ellipsoids found
} // end draw rand pixels in input ellipsoids

// draw 2d projections read from the JSON file at class github
function drawInputEllipsoidsUsingArcs(context) {
    var inputEllipsoids = getInputEllipsoids();
    
    
    if (inputEllipsoids != String.null) { 
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var w = context.canvas.width;
        var h = context.canvas.height;
        var n = inputEllipsoids.length; 
        //console.log("number of ellipsoids: " + n);

        // Loop over the ellipsoids, draw each in 2d
        for (var e=0; e<n; e++) {
            context.fillStyle = 
                "rgb(" + Math.floor(inputEllipsoids[e].diffuse[0]*255)
                +","+ Math.floor(inputEllipsoids[e].diffuse[1]*255)
                +","+ Math.floor(inputEllipsoids[e].diffuse[2]*255) +")"; // diffuse color
            context.save(); // remember previous (non-) scale
            context.scale(1, inputEllipsoids[e].b/inputEllipsoids[e].a); // scale by ellipsoid ratio 
            context.beginPath();
            context.arc(
                Math.round(w*inputEllipsoids[e].x),
                Math.round(h*inputEllipsoids[e].y),
                Math.round(w*inputEllipsoids[e].a),
                0,2*Math.PI);
            context.restore(); // undo scale before fill so stroke width unscaled
            context.fill();
            //console.log(context.fillStyle);
            //console.log("x: "+Math.round(w*inputEllipsoids[e].x));
            //console.log("y: "+Math.round(h*inputEllipsoids[e].y));
            //console.log("a: "+Math.round(w*inputEllipsoids[e].a));
            //console.log("b: "+Math.round(h*inputEllipsoids[e].b));
        } // end for ellipsoids
    } // end if ellipsoids found
} // end draw input ellipsoids

//put random points in the triangles from the class github
function drawRandPixelsInInputTriangles(context) {
    var inputTriangles = getInputTriangles();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputTriangles != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var cx = 0; var cy = 0; // init center x and y coord
        var numTrianglePixels = 0; // init num pixels in triangle
        var c = new Color(0,255,0,0); // init the triangle color
        var n = inputTriangles.length; // the number of input files
        //console.log("number of files: " + n);

        // Loop over the triangles, draw rand pixels in each
        for (var f=0; f<n; f++) {
        	var tn = inputTriangles[f].triangles.length;
        	//console.log("number of triangles in this files: " + tn);
        	
        	// Loop over the triangles, draw each in 2d
        	for(var t=0; t<tn; t++){
        		var vertex1 = inputTriangles[f].triangles[t][0];
        		var vertex2 = inputTriangles[f].triangles[t][1];
        		var vertex3 = inputTriangles[f].triangles[t][2];

        		var vertexPos1 = inputTriangles[f].vertices[vertex1];
        		var vertexPos2 = inputTriangles[f].vertices[vertex2];
        		var vertexPos3 = inputTriangles[f].vertices[vertex3];
        		//console.log("vertexPos1 " + vertexPos1);
        		//console.log("vertexPos2 " + vertexPos2);
        		//console.log("vertexPos3 " + vertexPos3);
        		
        		// triangle position on canvas
        		
        		var v1 = [w*vertexPos1[0], h*vertexPos1[1]];
        		var v2 = [w*vertexPos2[0], h*vertexPos2[1]];
        		var v3 = [w*vertexPos3[0], h*vertexPos3[1]];
        		
        		// calculate triangle area on canvas (shoelace formula)
        		var triangleArea = 0.5*Math.abs(v1[0]*v2[1]+v2[0]*v3[1]+v3[0]*v1[1]-v2[0]*v1[1]-v3[0]*v2[1]-v1[0]*v3[1]);
        		var numTrianglePixels = triangleArea; // init num pixels in triangle
            	//console.log("triangle area " + triangleArea);
            	numTrianglePixels *= PIXEL_DENSITY; // percentage of triangle area to render to pixels
            	numTrianglePixels = Math.round(numTrianglePixels);
            	// console.log("numTrianglePixels " + numTrianglePixels);
            	c.change(
            		inputTriangles[f].material.diffuse[0]*255,
                	inputTriangles[f].material.diffuse[1]*255,
                	inputTriangles[f].material.diffuse[2]*255,
                	255); // triangle diffuse color
            	for (var p=0; p<numTrianglePixels; p++) {
                    var point; // on canvas plane
            		var triangleTest = 0;
            		while (triangleTest == 0 ){ //if the pixel outside the triangle
                  
            			point = [Math.floor(Math.random()*w), Math.floor(Math.random()*h)];
                    	// plane checking
            			
                    	var t1 = ((point[0]-v2[0]) * (v1[1] - v2[1]) - (v1[0] - v2[0]) * (point[1] - v2[1])) < 0.0;
                    	var t2 = ((point[0]-v3[0]) * (v2[1] - v3[1]) - (v2[0] - v3[0]) * (point[1] - v3[1])) < 0.0;
                    	var t3 = ((point[0]-v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (point[1] - v1[1])) < 0.0;
                    	
                    	if((t1==t2)&&(t2==t3)) // draw the pixel if inside the triangle
                    		triangleTest = 1;
            		}
            		drawPixel(imagedata,point[0],point[1],c);
                	//console.log("color: ("+c.r+","+c.g+","+c.b+")");
                	//console.log("x: "+ x);
                	//console.log("y: "+ y);
            	} // end for pixels in triangle
        	} // end for triangles
    	} // end for files
        context.putImageData(imagedata, 0, 0);
    } // end if triangle file found
} // end draw rand pixels in input triangles

//draw 2d projections traingle from the JSON file at class github
function drawInputTrainglesUsingPaths(context) {
    var inputTriangles = getInputTriangles();
    
    if (inputTriangles != String.null) { 
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var w = context.canvas.width;
        var h = context.canvas.height;
        var n = inputTriangles.length; 
        //console.log("number of files: " + n);

        // Loop over the input files
        for (var f=0; f<n; f++) {
        	var tn = inputTriangles[f].triangles.length;
        	//console.log("number of triangles in this files: " + tn);
        	
        	// Loop over the triangles, draw each in 2d
        	for(var t=0; t<tn; t++){
        		var vertex1 = inputTriangles[f].triangles[t][0];
        		var vertex2 = inputTriangles[f].triangles[t][1];
        		var vertex3 = inputTriangles[f].triangles[t][2];

        		var vertexPos1 = inputTriangles[f].vertices[vertex1];
        		var vertexPos2 = inputTriangles[f].vertices[vertex2];
        		var vertexPos3 = inputTriangles[f].vertices[vertex3];
        		//console.log("vertexPos1 " + vertexPos1);
        		//console.log("vertexPos2 " + vertexPos2);
        		//console.log("vertexPos3 " + vertexPos3);
        		
            	context.fillStyle = 
            	    "rgb(" + Math.floor(inputTriangles[f].material.diffuse[0]*255)
            	    +","+ Math.floor(inputTriangles[f].material.diffuse[1]*255)
            	    +","+ Math.floor(inputTriangles[f].material.diffuse[2]*255) +")"; // diffuse color
            
            	var path=new Path2D();
            	path.moveTo(w*vertexPos1[0],h*vertexPos1[1]);
            	path.lineTo(w*vertexPos2[0],h*vertexPos2[1]);
            	path.lineTo(w*vertexPos3[0],h*vertexPos3[1]);
            	path.closePath();
            	context.fill(path);

        	} // end for triangles
        } // end for files
    } // end if triangle files found
} // end draw input triangles

// put random points in the boxes from the class github
function drawRandPixelsInInputBoxes(context) {
    var inputBoxes = getInputBoxes();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 100;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputBoxes != String.null) { 
	    var x  = 0; var y  = 0; // pixel coord init
        var lx = 0; var rx = 0; // input lx, rx from boxes.json
        var by = 0; var ty = 0; // input by, ty from boxes.json
        var fz = 0; var rz = 0; // input fz, rz from boxes.json
        var numBoxPixels = 0; // init num pixels in boxes
        var c = new Color(255,255,255,255); // init the box color
        var n = inputBoxes.length; // the number of input boxes
        //console.log("number of ellipses: " + n);

        // Loop over the ellipsoids, draw rand pixels in each
        for (var b=0; b<n; b++) {
			// input lx,rx,by,ty on canvas
			lx = w*inputBoxes[b].lx;
			rx = w*inputBoxes[b].rx;
			by = h*inputBoxes[b].by;
			ty = h*inputBoxes[b].ty;           
			
            numBoxesPixels  = (rx-lx)*(ty-by); // projected box area 
            numBoxesPixels *= PIXEL_DENSITY;  // percentage of box area to render to pixels
            numBoxesPixels  = Math.round(numBoxesPixels);
           
            //console.log("num box pixels: "+numBoxesPixels);
            
			c.change(
                inputBoxes[b].diffuse[0]*255,
                inputBoxes[b].diffuse[1]*255,
                inputBoxes[b].diffuse[2]*255,
                255); // box diffuse color
            for (var p=0; p<numBoxesPixels; p++) {
                do {
                    x = Math.floor(Math.random()*w); 
                    y = Math.floor(Math.random()*h); 
                } while ( x<lx || x>rx || y>ty || y<by ) // inside the projection
                drawPixel(imagedata,x,y,c);
                //console.log("color: ("+c.r+","+c.g+","+c.b+")");
                //console.log("x: " + x);
                //console.log("y: " + y);
            } // end for pixels in box
        } // end for boxes
        context.putImageData(imagedata, 0, 0);
    } // end if boxes found
} // end draw rand pixels in input boxes

//draw 2d projections boxes from the JSON file at class github
function drawInputBoxesUsingPaths(context) {
    var inputBoxes = getInputBoxes();
    var n = inputBoxes.length; // the number of input boxes
	
    if (inputBoxes != String.null) { 
		var w = context.canvas.width;
        var h = context.canvas.height;
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var x  = 0; var y  = 0; // pixel coord init
        var lx = 0; var rx = 0; // input lx, rx from boxes.json
        var by = 0; var ty = 0; // input by, ty from boxes.json
        var fz = 0; var rz = 0; // input fz, rz from boxes.json
        //console.log("number of files: " + n);

        // Loop over the input files
        for (var b=0; b<n; b++) {
				
			// input lx,rx,by,ty on canvas
			lx = w*inputBoxes[b].lx;
			rx = w*inputBoxes[b].rx;
			by = h*inputBoxes[b].by;
			ty = h*inputBoxes[b].ty; 
        		
            context.fillStyle = 
            	"rgb(" + Math.floor(inputBoxes[b].diffuse[0]*255)
            	+","+ Math.floor(inputBoxes[b].diffuse[1]*255)
            	+","+ Math.floor(inputBoxes[b].diffuse[2]*255) +")"; // diffuse color
            
            var path=new Path2D();
            path.moveTo(lx,ty);
            path.lineTo(lx,by);
            path.lineTo(rx,by);
			path.lineTo(rx,ty);
            path.closePath();
            context.fill(path);

        } // end for files
    } // end if box files found
} // end draw input boxes

function renderTriangles(context) {
    
    var inputTriangles = getInputTriangles();
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w, h);
    var depthBuffer = new Array(w * h).fill(Infinity);
    

    if (inputTriangles != String.null) {
        var eye = {x: 0.5, y: 0.5, z: -0.5};
        var lookAt = {x: 0, y: 0, z: 1};
        var up = {x: 0, y: 1, z: 0};

        for (var x = 0; x < w; x++) {
            for (var y = 0; y < h; y++) {
                var ray = pixelToRay(x, y, w, h, eye, lookAt, up);
                // console.log(ray)
                var closestIntersection = null;
                var closestTriangle = null;
                var closestTriangleData = null;
                // console.log(ray.direction)

                for (var f = 0; f < inputTriangles.length; f++) {
                    for (var t = 0; t < inputTriangles[f].triangles.length; t++) {
                        var triangle = getTriangleVertices(inputTriangles[f], t);
                        // triangle[0][0] = -triangle[0][0];
                        // triangle[1][0] = -triangle[1][0];
                        // triangle[2][0] = -triangle[2][0];
                        var intersection = rayTriangleIntersection(ray, triangle);
                        
                        if (intersection && (!closestIntersection || intersection.distance < closestIntersection.distance)) {
                            closestIntersection = intersection;
                            // console.log("closest intersection: " + JSON.stringify(closestIntersection));
                            closestTriangle = triangle;
                            closestTriangleData = inputTriangles[f];
                        }
                    }
                }
                // console.log("closest intersection: " + closestIntersection)
                if (closestIntersection) {
                    var depthIndex = y * w + x;
                    if (closestIntersection.distance < depthBuffer[depthIndex]) {
                        depthBuffer[depthIndex] = closestIntersection.distance;
                        var color = calculateLighting(closestIntersection, closestTriangle, closestTriangleData.material);
                        drawPixel(imagedata, x, y, color);
                    }
                }
            }
        }
        
        context.putImageData(imagedata, 0, 0);
    }
}

function pixelToRay(x, y, w, h, eye, lookAt, up) {
    // Calculate the view plane
    var viewPlaneDistance = 0.5;
    var aspectRatio = w / h;
    var viewPlaneHeight = 1;
    var viewPlaneWidth = aspectRatio * viewPlaneHeight;

    // Calculate view plane coordinates
    var u = (x / w - 0.5) * viewPlaneWidth;
    var v = (0.5 - y / h) * viewPlaneHeight;

    // Calculate view space basis vectors
    var forward = normalize(subtract(lookAt, eye));
    var right = normalize(crossProduct(forward, up));
    var upVector = crossProduct(right, forward);

    // Calculate ray direction
    var rayDirection = normalize({
        x: u * right.x + v * upVector.x + viewPlaneDistance * forward.x,
        y: u * right.y + v * upVector.y + viewPlaneDistance * forward.y,
        z: u * right.z + v * upVector.z + viewPlaneDistance * forward.z
    });

    return { origin: eye, direction: rayDirection };
}

function getTriangleVertices(triangleData, index) {
    var vertexIndices = triangleData.triangles[index];
    return vertexIndices.map(i => triangleData.vertices[i]);
}

function rayTriangleIntersection(ray, triangle) {
    const EPSILON = 0.0000001;
    var vertex0 = triangle[0];
    var vertex1 = triangle[1];
    var vertex2 = triangle[2];

    var edge1 = subtract(vertex1, vertex0);
    var edge2 = subtract(vertex2, vertex0);
    var h = crossProduct(ray.direction, edge2);
    var a = dotProduct(edge1, h);

    if (a > -EPSILON && a < EPSILON) return null;

    var f = 1.0 / a;
    var s = subtract(ray.origin, vertex0);
    var u = f * dotProduct(s, h);

    if (u < 0.0 || u > 1.0) return null;

    var q = crossProduct(s, edge1);
    var v = f * dotProduct(ray.direction, q);

    if (v < 0.0 || u + v > 1.0) return null;

    var t = f * dotProduct(edge2, q);
    // console.log("t:", t);

    if (t > EPSILON) {
        var scaledDirection = scaleVector(ray.direction, t);
        // console.log("scaledDirection:", scaledDirection);
        var point = add(ray.origin, scaledDirection);
        // console.log("point:", point);
        return {
            point: point,
            distance: t,
            u: u,
            v: v
        };
    }

    return null;
}

function scaleVector(v, s) {
    return { x: v.x * s, y: v.y * s, z: v.z * s };
}

function calculateLighting(intersection, triangle, material) {
    var lightPosition = { x: -3, y: 1, z: -0.5 };
    var lightColor = { r: 1, g: 1, b: 1 };
    var ambientIntensity = 0.1;
    var diffuseIntensity = 0.7;
    var specularIntensity = 0.2;
    var shininess = 10;

    // Calculate surface normal
    var edge1 = subtract(triangle[1], triangle[0]);
    var edge2 = subtract(triangle[2], triangle[0]);
    var normal = normalize(crossProduct(edge1, edge2));

    // Calculate light direction
    var lightDir = normalize(subtract(lightPosition, intersection.point));

    // Ambient component
    var ambient = scale(material.diffuse, ambientIntensity);
    

    // Diffuse component
    var diffuseFactor = Math.max(dotProduct(normal, lightDir), 0);
    var diffuse = scale(material.diffuse, diffuseFactor * diffuseIntensity);
    

    // Specular component
    var viewDir = normalize(subtract({ x: 0.5, y: 0.5, z: -0.5 }, intersection.point));
    var halfwayDir = normalize(add(lightDir, viewDir));
    var specularFactor = Math.pow(Math.max(dotProduct(normal, halfwayDir), 0), shininess);
    var specular = scale(lightColor, specularFactor * specularIntensity);

    // Combine components
    var finalColor = add(add(ambient, diffuse), specular);

    const finalColour = clampColor(finalColor);

    return new Color(
        Math.min(finalColour[0] * 255, 255),
        Math.min(finalColour[1] * 255, 255),
        Math.min(finalColour[2] * 255, 255),
        255
      );
}

function clampColor(color) {
    return [
      Math.max(0, Math.min(1, color[0])),
      Math.max(0, Math.min(1, color[1])),
      Math.max(0, Math.min(1, color[2]))
    ];
}

// Vector operations
function add(a, b) {
    const normalize = (obj) => {
      if (Array.isArray(obj)) return obj;
      if ('r' in obj) return [obj.r, obj.g, obj.b];
      if ('x' in obj) return [obj.x, obj.y, obj.z];
      if ('0' in obj) return [obj[0], obj[1], obj[2]];
      throw new Error("Unsupported input type for add function");
    };
  
    const arrA = normalize(a);
    const arrB = normalize(b);
  
    return [arrA[0] + arrB[0], arrA[1] + arrB[1], arrA[2] + arrB[2]];
  }
function subtract(a, b) {
    if (Array.isArray(a)) {
      a = { x: a[0], y: a[1], z: a[2] };
    }
    if (Array.isArray(b)) {
      b = { x: b[0], y: b[1], z: b[2] };
    }
    const epsilon = 1e-6; // small value to account for precision errors
    return { x: a.x - b.x + epsilon, y: a.y - b.y + epsilon, z: a.z - b.z + epsilon };
  }

  function scale(a, s) {
    if (Array.isArray(a)) {
      return [a[0] * s, a[1] * s, a[2] * s];
    } else if (typeof a === 'object' && 'r' in a) {
      return { r: a.r * s, g: a.g * s, b: a.b * s };
    } else {
      throw new Error("Unsupported input type for scale function");
    }
  }
function dotProduct(a, b) {
    if (Array.isArray(a)) {
      a = { x: a[0], y: a[1], z: a[2] };
    }
    if (Array.isArray(b)) {
      b = { x: b[0], y: b[1], z: b[2] };
    }
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

function crossProduct(a, b) {
    return {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x
    };
}

function normalize(v) {
    // console.log(v);
    if (Array.isArray(v)) {
        var length = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        // console.log("length: ", length);
        return [v[0] / length, v[1] / length, v[2] / length];
    } else {
        var length = Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        // console.log("Normalized object:", { x: v.x / length, y: v.y / length, z: v.z / length });
        return { x: v.x / length, y: v.y / length, z: v.z / length };
    }
}

/* main -- here is where execution begins after window load */

function main() {

    // Get the canvas and context
    var canvas = document.getElementById("viewport"); 
    var context = canvas.getContext("2d");
    renderTriangles(context);

    // drawRandPixelsInInputBoxes(context);
      // shows how to draw pixels and read input file
 
    // Create the image
    //drawRandPixels(context);
      // shows how to draw pixels
    
    //drawRandPixelsInInputEllipsoids(context);
      // shows how to draw pixels and read input file
      
    //drawInputEllipsoidsUsingArcs(context);
      // shows how to read input file, but not how to draw pixels
    
    // drawRandPixelsInInputTriangles(context);
      // shows how to draw pixels and read input file
    
    // drawInputTrainglesUsingPaths(context);
      // shows how to read input file, but not how to draw pixels
    
    //drawInputBoxesUsingPaths(context);
      // shows how to read input file, but not how to draw pixels
}
