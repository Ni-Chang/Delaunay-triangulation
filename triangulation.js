/**
 TODO Replace this by your own, correct, triangulation function
 Triangles should be return as arrays of array of indexes
 e.g., [[1,2,3],[2,3,4]] encodes two triangles, where the indices are relative to the array points
**/

var edges = [];
var faces = [];
// Left Right
var auxiliaryFace = [];
var auxiliaryEdge = -1;

/*   Find Super Triangle  */

function det(p, q, r){
	var detvalue = (q.x-p.x)*(r.y-p.y)-(r.x-p.x)*(q.y-p.y);
	if (detvalue > 0)
		return 1;
	if (detvalue == 0)
		return 0;
	return -1;
}

function det2d(u, v, w, z){
	var a1 = v.x-u.x;
	var b1 = w.x-u.x;
	var c1 = z.x-u.x;
	var a2 = v.y-u.y;
	var b2 = w.y-u.y;
	var c2 = z.y-u.y;
	var a3 = a1*(v.x+u.x) + a2*(v.y+u.y);
	var b3 = b1*(w.x+u.x) + b2*(w.y+u.y);
	var c3 = c1*(z.x+u.x) + c2*(z.y+u.y);
	var det2dresult = a1*b2*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1 - a1*b3*c2 - a2*b1*c3;
	if (det2dresult > 0)
		return 1;
	if (det2dresult < 0)
		return -1;
	return 0;

}

function findyMaxyMin(points){
	var yExtents = {"ymax": points[0].y, "ymaxX":points[0].x, "ymin":points[0].y}
	for (i=1; i<points.length; i++){
		if(points[i].y > yExtents.ymax){
			yExtents.ymax = points[i].y;
			yExtents.ymaxX = points[i].x;
		}
		if(points[i].y < yExtents.ymin){
			yExtents.ymin = points[i].y;
		}
	}
	yExtents.ymin--;
	return yExtents;
}

function findBigTriangle(points){
	yExtents = findyMaxyMin(points);
	var bigTriangle = [];
	for(i = 0; i<3; i++){
		bigTriangle.push(new Object())
	}
	bigTriangle[0] = {'x': yExtents.ymaxX-1, 'y': yExtents.ymax+1, 'z': 0.0};
	bigTriangle[1] = Object.assign({}, points[0]);
	bigTriangle[2] = Object.assign({}, points[0]);
	for(i = 1; i<points.length; i++){
		if(det(bigTriangle[0], bigTriangle[1], points[i]) < 0){
			bigTriangle[1] = Object.assign({}, points[i]);
		}
		if(det(bigTriangle[0], bigTriangle[2], points[i]) > 0){
			bigTriangle[2] = Object.assign({}, points[i]);
		}
	}
	bigTriangle[1].x = (yExtents.ymin - bigTriangle[1].y)/(bigTriangle[1].y-bigTriangle[0].y)*(bigTriangle[1].x-bigTriangle[0].x)+bigTriangle[1].x-1;
	bigTriangle[2].x = (yExtents.ymin - bigTriangle[2].y)/(bigTriangle[2].y-bigTriangle[0].y)*(bigTriangle[2].x-bigTriangle[0].x)+bigTriangle[2].x+1;
	bigTriangle[1].y = yExtents.ymin;
	bigTriangle[2].y = yExtents.ymin;
	points.push(bigTriangle[0],bigTriangle[1], bigTriangle[2]);
	faces.push(new face(0));
	faces.push(new face(0));
	edges.push(new edge(points.length-3, points.length-2, 1, 0, 2, 1));
	edges.push(new edge(points.length-2, points.length-1, 1, 0, 0, 2));
	edges.push(new edge(points.length-1, points.length-3, 1, 0, 1, 0));
	console.log(bigTriangle);
	return;
}

function ifccw(a, b, c){
	if(det(a, b, c) > 0)
		return 1;
	return 0;
}

/* Triangluation */
// 0-Not Intersect, 1-Intersect, -1-Interset at s1's endpoint, 2-shared point
function classifyIntersection(s1from, s1to, s2from, s2to) {
	var intersectionType, intersectionTypeDescription;
	// Compute the 4 determinal of the 2 segments
	det_s1f_s1t_s2f = det(s1from, s1to, s2from);
	det_s1f_s1t_s2t = det(s1from, s1to, s2to);
	det_s2f_s2t_s1f = det(s2from, s2to, s1from);
	det_s2f_s2t_s1t = det(s2from, s2to, s1to);
	determinal1 = det_s1f_s1t_s2f * det_s1f_s1t_s2t;
	determinal2 = det_s2f_s2t_s1f * det_s2f_s2t_s1t;

	// 2 points of s1 in the same side of s2, or 2 points of s2 in the same side of s1
	if (determinal1>0 || determinal2>0) {
		//Not Intersect
		intersectionType = 0;
	}
	// Both segments span the other one
	else if (determinal1 < 0 && determinal2 < 0) {
		//Intersect -- Intersection is an interior point
		intersectionType = 1;
	}
	// Complex situation
	else if (determinal1 == 0 && determinal2 == 0) {
			// Determine the positional relationship of the endpoints using x coordinates
		if (s1from.x != s2from.x) {
			s1_max = Math.max(s1from.x, s1to.x);
			s1_min = Math.min(s1from.x, s1to.x);
			s2_max = Math.max(s2from.x, s2to.x);
			s2_min = Math.min(s2from.x, s2to.x);
		}
		// If x coordinates are same, use y coordinates
		else {
			s1_max = Math.max(s1from.y, s1to.y);
			s1_min = Math.min(s1from.y, s1to.y);
			s2_max = Math.max(s2from.y, s2to.y);
			s2_min = Math.min(s2from.y, s2to.y);
		}
		// Beginning of s1 is in the left(top) of the beginning of s2
		if (s1_min < s2_min) {
			// End of s1 is in the left(top) of the beginning of s2
			if (s1_max < s2_min) {
				// Not Intersect
				intersectionType = 0;
			}
			else if (s1_max == s2_min) {
				intersectionType = 2; //shared endpoint
			}
			// End of s1 is on s2
			else if (s1_max <= s2_max){
				//Intersect -- Intersection is a segment;
				intersectionType = 1;
			}
		}
		// Beginning of s2 is in the left(top) of the beginning of s1
		else if (s1_min > s2_min) {
			// End of s2 is in the left(top) of the beginning of s1
			if(s2_max < s1_min){
				//Not Intersect
				intersectionType = 0;
			}
			else if(s2_max == s1_min) {
				intersectionType = 2;//shared endpoint
			}
			// End of s2 is on s1
			else if (s2_max <= s1_max){
				//Intersect -- Intersection is a segment
				intersectionType = 1;
			}
		}
	}
	// s2 has a endpoint on s1
	else if (determinal1 == 0) {
		//Intersect -- Intersection is s2's endpoint
		intersectionType = 1;
	}
	// s1 has a endpoint on s2
	else if (determinal2 == 0) {
		if(det_s2f_s2t_s1t==0){
			//query point is on boundary
			intersectionType = -1;
		}
		else{
			//Intersect -- Intersection is s1's endpoint
			intersectionType = 1;
		}
	}
	return intersectionType;
}

// 1-inside, -1-not inside, 0-boundary 2-vertex
function insideTriangle(p, triangle) {
	var det1 = det(triangle[0], triangle[1], p);
	var det2 = det(triangle[1], triangle[2], p);
	var det3 = det(triangle[2], triangle[0], p);
	var other2;
	// The point is on the same side of the 3 lines
	if ((det1 > 0 && det2 > 0 && det3 > 0) || (det1 < 0 && det2 < 0 && det3 < 0)) {
		return 1;//Inside
	}
	// The point is on a line or its extension
	else if((det1 * det2 * det3) == 0) {
		if (det1 == 0) {
			other2 = det2 * det3;
		}
		else if (det2 == 0) {
			other2 = det1 * det3;
		}
		else {
			other2 = det1 * det2;
		}
		// The point is on one line
		if (other2 > 0) {
			return 0;//Boundary
		}
		// The point is on two lines
		else if(other2 == 0){
			return 2;//Vertex
		}
		// The point is on the extension
		else {
			return -1;//Outside
		}
	}
	// The point is out of the triangle
	else{
		return -1;//Outside
	}
}

// 1-interior, 0-boundary or exterior
function insideCircle(p, circle_points) {
	// Let the three points on the circle be counterclockwise
	if(ifccw(circle_points[0], circle_points[1], circle_points[2]) == 1){
		var a = circle_points[0];
		var b = circle_points[1];
		var c = circle_points[2];
	}
	else{
		var a = circle_points[2];
		var b = circle_points[1];
		var c = circle_points[0];
	}
	
	var dirction = det2d(a, b, c, p);
	if (dirction < 0){
		return 1;//"Interior";
	}
	else if (dirction > 0){
		return 0;//"Exterior";
	}
	else{
		return 0;//"Boundary";
	}
}

function edge(vs, ve, fl, fr, ep, en){
	this.vs = vs;
	this.ve = ve;
	this.fl = fl;
	this.fr = fr;
	this.ep = ep;
	this.en = en;
}

function face(e){
	this.e = e;
}

function isCCW(edge, face){
	if(edges[edge].fl == face){
		return 1;
	}
	if(edges[edge].fr == face){
		return 0;
	}
	return -1;
}

function findTriEdges(face){
	var triEdges = new Array(3);
	triEdges[0] = faces[face].e;
	if(isCCW(triEdges[0], face) == 1){
		triEdges[1] = edges[triEdges[0]].ep;
	}
	else{
		triEdges[1] = edges[triEdges[0]].en;
	}
	if(isCCW(triEdges[1], face) == 1){
		triEdges[2] = edges[triEdges[1]].ep;
	}
	else{
		triEdges[2] = edges[triEdges[1]].en;
	}
	return triEdges;
}

function findTriPoints(refface){
	var triEdges = new Array(3);
	var triPoints = new Array(3);

	triEdges[0] = faces[refface].e;
	if(isCCW(triEdges[0], refface) == 1){
		triPoints[0] = edges[triEdges[0]].ve;
		triPoints[1] = edges[triEdges[0]].vs;
		triEdges[1] = edges[triEdges[0]].ep;
	}
	else{
		triPoints[0] = edges[triEdges[0]].vs;
		triPoints[1] = edges[triEdges[0]].ve;
		triEdges[1] = edges[triEdges[0]].en;
	}
	if(isCCW(triEdges[1], refface) == 1){
		triPoints[2] = edges[triEdges[1]].vs;
	}
	else{
		triPoints[2] = edges[triEdges[1]].ve;
	}
	return triPoints;
}

// CW from face's edge
// Add 2 faces and 3 edges
// update 3 edges
function updateFacesAndEdges(refface, point){
	var triEdges = new Array(3);
	var triPoints = new Array(3);

	triEdges[0] = faces[refface].e;
	if(isCCW(triEdges[0], refface) == 1){
		triPoints[0] = edges[triEdges[0]].ve;
		triPoints[1] = edges[triEdges[0]].vs;
		triEdges[1] = edges[triEdges[0]].ep;
	}
	else{
		triPoints[0] = edges[triEdges[0]].vs;
		triPoints[1] = edges[triEdges[0]].ve;
		triEdges[1] = edges[triEdges[0]].en;
	}
	if(isCCW(triEdges[1], refface) == 1){
		triEdges[2] = edges[triEdges[1]].ep;
		triPoints[2] = edges[triEdges[1]].vs;
	}
	else{
		triEdges[2] = edges[triEdges[1]].en;
		triPoints[2] = edges[triEdges[1]].ve;
	}
	
	faces.push(new face(triEdges[1]));
	faces.push(new face(triEdges[2]));
	edges.push(new edge(point, triPoints[0], faces.length-1, refface, edges.length+2, triEdges[0]));
	edges.push(new edge(point, triPoints[1], refface, faces.length-2, edges.length-1, triEdges[1]));
	edges.push(new edge(point, triPoints[2], faces.length-2, faces.length-1, edges.length-1, triEdges[2]));
	if(isCCW(triEdges[0], refface) == 1){
		edges[triEdges[0]].ep = edges.length-2;
	}
	else{
		edges[triEdges[0]].en = edges.length-2;
	}
	if(isCCW(triEdges[1], refface) == 1){
		edges[triEdges[1]].ep = edges.length-1;
		edges[triEdges[1]].fl = faces.length-2;
	}
	else{
		edges[triEdges[1]].en = edges.length-1;
		edges[triEdges[1]].fr = faces.length-2;
	}
	if(isCCW(triEdges[2], refface) == 1){
		edges[triEdges[2]].ep = edges.length-3;
		edges[triEdges[2]].fl = faces.length-1;
	}
	else{
		edges[triEdges[2]].en = edges.length-3;
		edges[triEdges[2]].fr = faces.length-1;
	}
	return;
}

function updateAuxiliary(refTriPoints, point){
	var insideTri1 = insideTriangle(points[0], [points[refTriPoints[0]], points[refTriPoints[2]], point]);
	var insideTri2 = insideTriangle(points[0], [points[refTriPoints[1]], points[refTriPoints[2]], point]);
	if(insideTri1 == 1){
		auxiliaryFace[0] = faces.length-1;
	}
	else if(insideTri2 == 1){
		auxiliaryFace[0] = faces.length-2;
	}
	if(auxiliaryEdge == -1){
		if(insideTri1 == 0 || insideTri2 == 0){
			if(insideTri1 == 0 && insideTri2 == 0){
				auxiliaryEdge = edges.length-1;
			}
			else if(insideTri1 == 0){
				auxiliaryEdge = edges.length-3;
			}
			else if(insideTri2 == 0){
				auxiliaryEdge = edges.length-2;
			}
			auxiliaryFace[0] = edges[auxiliaryEdge].fl;
			auxiliaryFace[1] = edges[auxiliaryEdge].fr;
		}
	}
	return;
}

function pointProcessing(pointIndex, findFace, enteredEdge){
	var point = points[pointIndex];
	var k = -1;
	var triEdges = findTriEdges(findFace);
	for(var i=0; i<3; i++){
		if(triEdges[i]!=enteredEdge){
			// 0-Not Intersect, 1-Intersect, -1-Interset at s1's endpoint
			var Intersection = classifyIntersection(points[0], point, points[edges[triEdges[i]].vs], points[edges[triEdges[i]].ve]);
			if(Intersection == 2){
				return [];
			}
			if(Intersection == 1){
				k = i;
				break;
			}
			else if(Intersection == -1){
				k = i + 3;
			}
		}
	}
	var refTriPoints = findTriPoints(findFace);
	if(k == -1){
		updateFacesAndEdges(findFace, pointIndex);
		updateAuxiliary(refTriPoints, point);
		return triEdges;
	}
	else if(k > 2){
		var coincideEdge = refTriPoints[k-3];
		var twoTriEdges = new Array(2);
		twoTriEdges[0] = findTriEdges(edges[coincideEdge].fl);
		twoTriEdges[1] = findTriEdges(edges[coincideEdge].fr);
		var fourEdges = new Array();
		for(var i = 0 ; i < 2; i++){
			for(var j = 0; j < 3; j++){
				if(twoTriEdges[i][j]!=coincideEdge){
					fourEdges.push(twoTriEdges[i][j])
				}
			}
		}
		var refTriPointsLeft = findTriPoints(edges[coincideEdge].fl);
		var refTriPointsRight = findTriPoints(edges[coincideEdge].fr);
		var thirdPointLeft;
		var thirdPointRight;
		for(var i=0; i<3; i++){
			if(refTriPointsLeft[i]!=edges[coincideEdge].vs && refTriPointsLeft[i]!=edges[coincideEdge].ve){
				thirdPointLeft = refTriPointsLeft[i];
			}
			if(refTriPointsRight[i]!=edges[coincideEdge].vs && refTriPointsRight[i]!=edges[coincideEdge].ve){
				thirdPointRight = refTriPointsRight[i];
			}
		}
		var leftThreeEdges = findTriEdges(edges[coincideEdge].fl);
		var rightThreeEdges = findTriEdges(edges[coincideEdge].fr);
		for(var i=0; i<3; i++){
			if(edges[leftThreeEdges[i]].ep == coincideEdge || edges[leftThreeEdges[i]].en == coincideEdge){
				var tempnext1 = leftThreeEdges[i];
			}
			if(edges[rightThreeEdges[i]].ep == coincideEdge || edges[rightThreeEdges[i]].en == coincideEdge){
				var tempnext2 = rightThreeEdges[i];
			}
		}
		faces.push(new face(edges.length+1));
		faces.push(new face(edges.length+2));
		edges.push(new edge(edges[coincideEdge].vs, pointIndex, faces.length-2, faces.length-1, edges[coincideEdge].ep, edges.length+2));
		edges[coincideEdge].vs = pointIndex;
		edges[coincideEdge].ep = edges.length;
		edges.push(new edge(pointIndex, thirdPointLeft, faces.length-2, edges[coincideEdge].fl, edges.length-1, tempnext1));
		edges.push(new edge(pointIndex, thirdPointRight, edges[coincideEdge].fr, faces.length-1, coincideEdge, tempnext2));
		
		if(edges[edges[edges.length-3].ep].ep == tempnext1){
			edges[edges[edges.length-3].ep].ep = edges.length-2;
			edges[edges[edges.length-3].ep].fl = faces.length-2;
		}
		else{
			edges[edges[edges.length-3].ep].en = edges.length-2;
			edges[edges[edges.length-3].ep].fr = faces.length-2;
		}
		if(edges[edges[coincideEdge].en].ep == tempnext2){
			edges[edges[coincideEdge].en].ep = edges.length-1;
		}
		else{
			edges[edges[coincideEdge].en].en = edges.length-1;
		}
		if(edges[tempnext2].ep == coincideEdge){
			edges[tempnext2].fl = faces.length-1;
		}
		else{
			edges[tempnext2].fr = faces.length-1;
		}
		//(coincideFaces[j], pointIndex);
		updateAuxAfterFlip(edges[coincideEdge].fl, faces.length-2);
		updateAuxAfterFlip(edges[coincideEdge].fr, faces.length-1);
		//coincideFaces[j].e = -1;
		return fourEdges;
	}
	else{
		var newFace;
		if(edges[triEdges[k]].fl == findFace){
			newFace = edges[triEdges[k]].fr;
		}
		else{
			newFace = edges[triEdges[k]].fl;
		}
		return pointProcessing(pointIndex, newFace, triEdges[k]);
	}
}

function findPointNotOnEdge(faceIndex, edgeIndex){
	var triPoints = findTriPoints(faceIndex);
	var vs = edges[edgeIndex].vs;
	var ve = edges[edgeIndex].ve;
	for(var i=0; i<3; i++){
		if(triPoints[i]!=vs && triPoints[i]!=ve){
			return triPoints[i];
		}
	}
}

function updateAuxAfterFlip(face1, face2, flipedEdge){
	var triPointsIndex1 = findTriPoints(face1);
	var triPointsIndex2 = findTriPoints(face2);
	var triPoints1 = new Array(3);
	var triPoints2 = new Array(3);
	for(var i=0; i<3; i++){
		triPoints1[i] = points[triPointsIndex1[i]];
		triPoints2[i] = points[triPointsIndex2[i]];
	}
	// 1-inside, -1-not inside, 0-boundary 2-vertex
	if(auxiliaryEdge==-1){
		if(insideTriangle(points[0], triPoints1) == 1){
			auxiliaryFace[0] = face1;
		}
		else if(insideTriangle(points[0], triPoints2) == 1){
			auxiliaryFace[0] = face2;
		}
		else if(insideTriangle(points[0], triPoints1)==2 || insideTriangle(points[0], triPoints2)==2){
			auxiliaryEdge = flipedEdge;
			auxiliaryFace[0] = edges[flipedEdge].fl;
			auxiliaryFace[1] = edges[flipedEdge].fr;
		}
	}
	else{
		if(insideTriangle(points[0], triPoints1) == 1){
			auxiliaryFace[0] = face1;
			auxiliaryEdge=-1;
		}
		else if(insideTriangle(points[0], triPoints2) == 1){
			auxiliaryFace[0] = face2;
			auxiliaryEdge=-1;
		}
	}
	return;
}

function delaunayFlip(pointIndex, checkEdge){
	var outsideFace;
	if(det(points[edges[checkEdge].vs], points[edges[checkEdge].ve], points[pointIndex])>0){
		outsideFace = edges[checkEdge].fr;
	}
	else{
		outsideFace = edges[checkEdge].fl;
	}
	if(outsideFace!=0){
		outsideTriPointsIndex = findTriPoints(outsideFace);
		var outsideTriPoints = new Array(3);
		for (var i=0; i<3; i++){
			outsideTriPoints[i] = points[outsideTriPointsIndex[i]];
		}
		var flag = insideCircle(points[pointIndex], outsideTriPoints);
		if(flag == 1){
			faces[edges[checkEdge].fl].e = checkEdge;
			faces[edges[checkEdge].fr].e = checkEdge;
			var tempvs = findPointNotOnEdge(edges[checkEdge].fl, checkEdge);
			var tempve = findPointNotOnEdge(edges[checkEdge].fr, checkEdge);
			edges[checkEdge].vs = tempvs;
			edges[checkEdge].ve = tempve;
			var checkEdgePre = edges[checkEdge].ep;
			var checkEdgeLeftTriEdges = findTriEdges(edges[checkEdge].fl);
			var checkEdgeNext = edges[checkEdge].en;
			var checkEdgeRightTriEdges = findTriEdges(edges[checkEdge].fr);
			for(var i=0; i<3; i++){
				if(checkEdgeLeftTriEdges[i]!=checkEdge && checkEdgeLeftTriEdges[i]!=checkEdgePre){
					edges[checkEdge].ep = checkEdgeLeftTriEdges[i];
					break;
				}
			}
			for(var i=0; i<3; i++){
				if(checkEdgeRightTriEdges[i]!=checkEdge && checkEdgeRightTriEdges[i]!=checkEdgeNext){
					edges[checkEdge].en = checkEdgeRightTriEdges[i];
					break;
				}
			}
			if(edges[checkEdgePre].ep==edges[checkEdge].ep){

				edges[checkEdgePre].ep = checkEdge;
				edges[checkEdgePre].fl = edges[checkEdge].fr;
			}
			else{
				edges[checkEdgePre].en = checkEdge;
				edges[checkEdgePre].fr = edges[checkEdge].fr;
			}
			if(edges[checkEdgeNext].ep==edges[checkEdge].en){
				edges[checkEdgeNext].ep = checkEdge;
				edges[checkEdgeNext].fl = edges[checkEdge].fl;
			}
			else{
				edges[checkEdgeNext].en = checkEdge;
				edges[checkEdgeNext].fr = edges[checkEdge].fl;
			}
			if(edges[edges[checkEdge].ep].ep == checkEdge){
				edges[edges[checkEdge].ep].ep = checkEdgeNext;
			}
			else{
				edges[edges[checkEdge].ep].en = checkEdgeNext;
			}
			if(edges[edges[checkEdge].en].ep == checkEdge){
				edges[edges[checkEdge].en].ep = checkEdgePre;
			}
			else{
				edges[edges[checkEdge].en].en = checkEdgePre;
			}
			var checkEdgeNewPre = edges[checkEdge].ep;
			var checkEdgeNewNext = edges[checkEdge].en;
			if(edges[checkEdge].vs == pointIndex){
				delaunayFlip(pointIndex, checkEdgeNext);
				delaunayFlip(pointIndex, checkEdgeNewNext);
			}
			else{
				delaunayFlip(pointIndex, checkEdgePre);
				delaunayFlip(pointIndex, checkEdgeNewPre);
			}
			updateAuxAfterFlip(edges[checkEdge].fl, edges[checkEdge].fr, checkEdge);
		}
	}
	return;
}

function computeTriangulation(points) {
	// Note that this does NOT return a triangulation, just a subset of it
	//var newPoints = points.sort(function(a,b) { if ((a.x) < (b.x)) return -1; else return 1;})
	findBigTriangle(points);
	auxiliaryFace[0] = 1;
	var outputTriangles = new Array(); 
	var triEdgesNow;
	for (var i=1; i<points.length-3; i++){
		if(auxiliaryEdge!=-1){
			if(det(edges[auxiliaryEdge].vs, edges[auxiliaryEdge].ve, points[i])>0){
				triEdgesNow = pointProcessing(i, auxiliaryFace[0], -1);
				for(var j=0; j<triEdgesNow.length; j++){
					delaunayFlip(i, triEdgesNow[j]);
				}
			}
			else{
				triEdgesNow = pointProcessing(i, auxiliaryFace[1], -1);
				for(var j=0; j<triEdgesNow.length; j++){
					delaunayFlip(i, triEdgesNow[j]);
				}
			}
		}
		else{
			triEdgesNow = pointProcessing(i, auxiliaryFace[0], -1);
			for(var j=0; j<triEdgesNow.length; j++){
				delaunayFlip(i, triEdgesNow[j]);
			}
		}	
		//console.log("point: "+i+"finished");
	}
	triEdgesNow = pointProcessing(0, auxiliaryFace[0], -1);
	for(var j=0; j<triEdgesNow.length; j++){
		delaunayFlip(0, triEdgesNow[j]);
	}

	//Output
	for (i=1; i<faces.length; i++){
		if(faces[i].e != -1){
			var triPoint1 = edges[faces[i].e].vs;
			var triPoint2 = edges[faces[i].e].ve;
			var triPoint3;
			if (edges[faces[i].e].fl == i){
				if (edges[edges[faces[i].e].ep].vs == triPoint1){
					triPoint3 = edges[edges[faces[i].e].ep].ve;
				}
				else{
					triPoint3 = edges[edges[faces[i].e].ep].vs;
				}
							
			}
			else{
				if (edges[edges[faces[i].e].en].vs == triPoint2){
					triPoint3 = edges[edges[faces[i].e].en].ve;
				}
				else{
					triPoint3 = edges[edges[faces[i].e].en].vs;
				}
					
			}
			if(det(points[triPoint1], points[triPoint2], points[triPoint3])>0){
				outputTriangles.push([triPoint1, triPoint2, triPoint3]);
			}
			else{
				outputTriangles.push([triPoint3, triPoint2, triPoint1]);
			}
		}
	}
	return outputTriangles;
}

function checkOutsidePolygon(p1Index, p2Index, p3Index, checkPointIndex){
	p1 = points[p1Index];
	p2 = points[p2Index];
	p3 = points[p3Index];
	checkPoint = points[checkPointIndex];
	if(det(p1, p2, p3) > 0){
		if(det(p1, p2, checkPoint)<0 || det(p2, p3, checkPoint)<0){
			return 1;
		}
		return 0;
	}
	else{
		if(det(p1, p2, checkPoint)<0 && det(p2, p3, checkPoint)<0){
			return 1;
		}
		return 0;
	}	
}

function pruneBoundaries(boundaries){
	for(var i=0; i<boundaries.length; i++){
		bdysize = boundaries[i].length;
		if(bdysize == 1){
			for(var k = 0; k<edges.length; k++){
				if(edges[k].vs == boundaries[i][0] || edges[k].ve == boundaries[i][0]){
					faces[edges[k].fl].e = -1;
					faces[edges[k].fr].e = -1;
					console.log(edges[k].fl+"deleted");
				}
			}
		}
		else{
			for(var k = 0; k<edges.length; k++){
				if(edges[k].vs == boundaries[i][0] || edges[k].ve == boundaries[i][0]){
					var checkPoint = (edges[k].vs == boundaries[i][0])?edges[k].ve:edges[k].vs;
					if(checkOutsidePolygon(boundaries[i][boundaries[i].length-1], boundaries[i][0], boundaries[i][1], checkPoint)){
						faces[edges[k].fl].e = -1;
						faces[edges[k].fr].e = -1;
						console.log(edges[k].fl+"deleted");
					}
				}
			}
			for(var j=1; j<boundaries[i].length-1; j++){
				for(var k = 0; k<edges.length; k++){
					if(edges[k].vs == boundaries[i][j] || edges[k].ve == boundaries[i][j]){
						var checkPoint = (edges[k].vs == boundaries[i][j])?edges[k].ve:edges[k].vs;
						if(checkOutsidePolygon(boundaries[i][j-1], boundaries[i][j], boundaries[i][j+1], checkPoint)){ 
							faces[edges[k].fl].e = -1;
							faces[edges[k].fr].e = -1;
							console.log(edges[k].fl+"deleted");
						}
					}
				}
			}
			for(var k = 0; k<edges.length; k++){
				if(edges[k].vs == boundaries[i][j] || edges[k].ve == boundaries[i][j]){
					var checkPoint = (edges[k].vs == boundaries[i][j])?edges[k].ve:edges[k].vs;
					if(checkOutsidePolygon(boundaries[i][j-1], boundaries[i][j], boundaries[i][0], checkPoint)){
						faces[edges[k].fl].e = -1;
						faces[edges[k].fr].e = -1;
						console.log(edges[k].fl+"deleted");
					}
				}
			}
		}
	}
	points.pop();
	points.pop();
	points.pop();
	var outputTriangles = new Array(); 
	for (i=1; i<faces.length; i++){
		if(faces[i].e != -1){
			var triPoint1 = edges[faces[i].e].vs;
			var triPoint2 = edges[faces[i].e].ve;
			var triPoint3;
			if (edges[faces[i].e].fl == i){
				if (edges[edges[faces[i].e].ep].vs == triPoint1){
					triPoint3 = edges[edges[faces[i].e].ep].ve;
				}
				else{
					triPoint3 = edges[edges[faces[i].e].ep].vs;
				}
							
			}
			else{
				if (edges[edges[faces[i].e].en].vs == triPoint2){
					triPoint3 = edges[edges[faces[i].e].en].ve;
				}
				else{
					triPoint3 = edges[edges[faces[i].e].en].vs;
				}
					
			}
			if(det(points[triPoint1], points[triPoint2], points[triPoint3])>0){
				outputTriangles.push([triPoint1, triPoint2, triPoint3]);
			}
			else{
				outputTriangles.push([triPoint3, triPoint2, triPoint1]);
			}
		}
	}
	return outputTriangles;
}