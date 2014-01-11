
/*

TODO:

	Move to using typed arrays for the particle properties

	Allow user to "draw" fixed particles on the field
		The particle classes would be created by the user on the page or from presets

	Importing/exporting of settings
		This could include the arrangement of "drawn" particles

	Concurrency?

*/

// To keep the compiler from complaining

/** @type {Object} */
var star;

/** @type {Object} */
var $;

(function() {
'use strict';

//
//  Body
//________//

// Represents a particle, with a position, radius, mass, acceleration, velocity, force, and color
// The color, size, and mass come from the type of particle

function Body(x, y, type) {
	// Set the position and type
	this.x = x;
	this.y = y;
	this.setType(type);
}

Body.prototype = {

	// The type or "class" of body
	type: null,

	// Whether it can move or not
	canMove: true,

	// The cell it's in
	cell: null,

	// Position, velocity, force, acceleration
	x: 0,
	y: 0,
	vx: 0,
	vy: 0,
	fx: 0,
	fy: 0,
	ax: 0,
	ay: 0,

	// Updates the kinematic properties of the body in regards to the specified amount of time
	applyForce: function Body_applyForce(dt) {
		this.ax = this.fx / this.mass;
		this.ay = this.fy / this.mass;
		this.vx += this.ax * dt;
		this.vy += this.ay * dt;
		this.x += this.vx * dt;
		this.y += this.vy * dt;
		// Drop the forces
		// ???: Why is this needed?
		this.fx = 0;
		this.fy = 0;
	},

	// Sets the body type
	setType: function Body_setType(type) {
		this.type = type;
		this.radius = star.classes[type].diameter / 2;
		this.mass = star.classes[type].mass;
		this.threshold = star.classes[type].threshold;
	}

};

//
//  Cell
//________//

// Represents a cell in the grid, to hold an arbitrary number of bodies

function Cell() {
	// Initiate the bodies array
	this.bodies = [];
}

Cell.prototype = {

	// Array of bodies
	bodies: null,

	// Origin of cell
	x: 0,
	y: 0,

	// Total mass in cell
	mass: 0,

	// Position of center of mass of cell
	massX: 0,
	massY: 0,

	// These were previously used as pointers to adjacent cells
	// Easy links to adjacent cells for cross-cell gravity and collisions
	// northernCell: null,
	// southernCell: null,
	// easternCell: null,
	// westernCell: null,

	// Clears the cell and sets it's position
	reset: function Cell_reset(x, y) {
		// Set the properties
		this.x = x;
		this.y = y;
		this.mass = 0;
		this.massX = 0;
		this.massY = 0;
		// Just clear the bodies array rather than create a new one
		this.bodies.length = 0;
	}

};

//
//  Star
//________//

window.star = {

	SIMULATION_INTERVAL_TIME: 20,
	//#replace star.SIMULATION_INTERVAL_TIME 20

	MAX_DISPLAY_WIDTH: 1920,
	//#replace star.MAX_DISPLAY_WIDTH 1920

	MAX_DISPLAY_HEIGHT: 1080,
	//#replace star.MAX_DISPLAY_HEIGHT 1080

	// The amount of time the physics simulation is allowed to account for
	// Lower values will make it progress slower but more accurately
	MAX_SIMULATION_STEP_TIME: 30,
	//#replace star.MAX_SIMULATION_STEP_TIME 30

	INIT_GRID_CELL_SIZE: 12000,
	//#replace star.INIT_GRID_CELL_SIZE 12000

	LEFT_CLICK_BODY_CLASS: 8,
	//#replace star.LEFT_CLICK_BODY_CLASS 8

		// Main

	particleCount: 3000,

	frameInterval: null,

	simTime: null,

	drawTime: null,

	enableAutomaticGridSizing: true,

	running: false,

	flipNextCollisionCheckDirection: false,

	bodies: [],

		// Physics

	gravity: 0.0025,

	wallBounce: -0.8,
	// Negative to account for bounce = -bounce * xVel

	wallFrict: 0.8,

	particleBounce: 0.7,

		// Rendering

	drawGrid: false,

	drawCellMass: false,

	drawCellCenterOfMass: false,

	particleDrawSizeMultiplier: 0.8,

		// Particle Class Data

	classes: [
		{},
		{
			color: 'Black',
			diameter: 3,
			mass: 20,
			threshold: 1000000,
			count: 0
		},
		{
			color: 'Purple',
			diameter: 4,
			mass: 30,
			threshold: 3*1000000,
			count: 0
		},
		{
			color: 'Blue',
			diameter: 5,
			mass: 40,
			threshold: 8*1000000,
			count: 0
		},
		{
			color: 'Green',
			diameter: 6,
			mass: 50,
			threshold: 12*1000000,
			count: 0
		},
		{
			color: 'Yellow',
			diameter: 7,
			mass: 60,
			threshold: 25*1000000,
			count: 0
		},
		{
			color: 'Orange',
			diameter: 8,
			mass: 70,
			threshold: 36*1000000,
			count: 0
		},
		{
			color: 'Red',
			diameter: 10,
			mass: 80,
			threshold: 49*1000000,
			count: 0
		},
		{
			color: 'Pink',
			diameter: 20,
			mass: -10000
		}
	],

	init: function star_init() {

		star.canvas = document.getElementById('canvas');

		star.setGridCellSize(star.INIT_GRID_CELL_SIZE);

		document.body.onresize = star.updateCanvasSize;
		star.updateCanvasSize();

		// Get the context
		star.ctx = star.canvas.getContext('2d');

		// Set a click handler to add a negative particle
		star.canvas.onmousedown = function(event) {
			var borderWidth = (star.canvas.offsetWidth - star.width)/2,
				b = new Body(
					event.clientX - star.canvas.offsetLeft - borderWidth,
					event.clientY - star.canvas.offsetTop - borderWidth,
					star.LEFT_CLICK_BODY_CLASS
				);
			b.canMove = false;
			star.bodies.push(b);
		};

		star.canvas.onmouseup = function() {
			star.bodies.pop();
		};

		star.ui.init();

		star.restart();

	},

	updateCanvasSize: function star_updateCanvasSize() {
		var container = document.getElementById('canvasContainer'),
			width = container.offsetWidth-40,
			height = container.offsetHeight-40;
		if(width > window.innerWidth) {
			width = window.innerWidth-40;
		}
		if(height > window.innerHeight) {
			height = window.innerHeight-40;
		}
		// Cap the width and height
		if(width > star.MAX_DISPLAY_WIDTH){
			width = star.MAX_DISPLAY_WIDTH;
		}
		if(height > star.MAX_DISPLAY_HEIGHT){
			height = star.MAX_DISPLAY_HEIGHT;
		}
		if(height < document.body.clientHeight-40) {
			height = document.body.clientHeight-40;
		}
		// Set the size
		star.width = width;
		star.height = height;
		star.canvas.width = width;
		star.canvas.height = height;
		// Reset the grid
		// star.setGridSize(star.columnCount, star.rowCount);
	},

	stepSimulation: function star_stepSimulation() {

		// Get the time elapsed since the last frame
		var now = Date.now(),
			elapsed = now-star.lastTime,
			s0;
		star.lastTime = now;

		if(star.running) {

			// Cap the simulation to keep it from being too inaccurate
			if(elapsed > star.MAX_SIMULATION_STEP_TIME) {
				elapsed = star.MAX_SIMULATION_STEP_TIME;
			}

			// Step the simulation
			s0 = Date.now();
			// Constrain the bodies
			star.moveConstrainsAll();
			// Apply the force on the bodies
			star.moveForceAll(elapsed);
			// Save the time taken
			star.simTime = Date.now()-s0;

			// Push the cell size and simulation time
			star.previousGridCellSizes.push(star.gridCellSize);
			star.previousSimTimes.push(star.simTime);
			if(star.previousGridCellSizes.length > star.PREVIOUS_GRID_CELL_SIZES_MAX_SIZE) {
				star.previousGridCellSizes.shift();
				star.previousSimTimes.shift();
			}

			// Render the scene
			s0 = Date.now();
			star.render();
			star.drawTime = Date.now()-s0;

		}

	},

	start: function star_start() {
		// Save the current time and start the frame interval
		star.lastTime = Date.now();
		star.frameInterval = setInterval(star.stepSimulation, star.SIMULATION_INTERVAL_TIME);
		// Setup the grid resizer
		if(star.enableAutomaticGridSizing) {
			star.startGridSizingRoutine();
		}
		star.running = true;
	},

	stop: function star_stop() {
		// Clear the frame interval and grid resizer
		clearInterval(star.frameInterval);
		star.clearGridSizingRoutine();
		star.running = false;
	},

	restart: function star_restart() {
		var i;
		// Set the default cell size
		star.setGridCellSize(star.defaultCellSize);
		// Clear the particle class counts
		for(i=0;i<star.classes.length;i++) {
			star.classes[i].count = 0;
		}
		star.classes[1].count = star.particleCount;
		// Fill with bodies
		star.bodies = [];
		for(i=star.particleCount-1; i>=0; i--) {
			star.bodies.push( new Body(Math.random()*star.width, Math.random()*star.height, 1) );
		}
		star.ui.updateClassCounts();
		star.start();

	},

	//
	//  Grid
	//________//

	cells: [],

	cellCount: 0,

	columns: null,

	rows: null,

	cellWidth: null,

	cellHeight: null,

	// needsResize: false,

	// newColumnCount: null,

	// newRowCount: null,

	// Inserts every body into the correct cell
	// Assumes that each particle cannot span more than one cell

	constrainAndSortBodies: function star_constrainAndSortBodies() {

		var cells = star.cells,
			bodies = star.bodies,
			i, body,
			x, y, index, cell,
			width = star.width,
			height = star.height,
			columns = star.columnCount,
			rows = star.rowCount,
			cellWidth = star.cellWidth,
			cellHeight = star.cellHeight,
			bounce = star.wallBounce,
			frict = star.wallFrict;

		// Resize if pending
		// if(star.needsResize) {
		// 	star.setSize(star.newColumnCount, star.newRowCount);
		// 	star.needsResize = false;
		// }

		// Reset all the cells
		for( i=star.cellCount-1; i>=0; i-- ) {
			cell = cells[i];
			cell.bodies.length = 0;
			cell.mass = 0;
			cell.massX = 0;
			cell.massY = 0;
		}

		// Loop over every body
		for( i=bodies.length-1; i>=0; i-- ) {
			body = bodies[i];

			// Constrain the body to within the field
			if(body.x < 0) {
				body.x = 0;
				body.vx = bounce*body.vx;
				body.vy *= frict;
			}
			else if(body.x > width) {
				body.x = width;
				body.vx = bounce*body.vx;
				body.vy *= frict;
			}

			if(body.y < 0) {
				body.y = 0;
				body.vy = bounce*body.vy;
				body.vx *= frict;
			}
			else if(body.y > height) {
				body.y = height;
				body.vy = bounce*body.vy;
				body.vx *= frict;
			}

			// Get the index of the cell it's in
			// Bitwise |0 is used to floor the number
			x = (body.x/cellWidth)|0;
			y = (body.y/cellHeight)|0;

			// In the event of x/y being at the absolute last row/column, knock it down by 1
			if(x === columns) {
				x--;
			}
			if(y === rows) {
				y--;
			}

			index = y*columns + x;
			cell = cells[index];
			cell.bodies.push(body);
			cell.mass += body.mass;
			cell.massX += body.x*body.mass;
			cell.massY += body.y*body.mass;
			body.cell = cell;

		}

		// Normalize each cell's center of mass 
		for( i=star.cellCount-1; i>=0; i-- ) {
			cell = cells[i];
			cell.massX = cell.massX / cell.mass;
			cell.massY = cell.massY / cell.mass;
		}

	},

	// pendResize: function star_pendResize(columns, rows) {
	// 	star.needsResize = true;
	// 	star.newColumnCount = columns;
	// 	star.newRowCount = rows;
	// },

	setGridSize: function star_setGridSize(columns, rows) {

		var y, x, cell;

		star.columnCount = columns;
		star.rowCount = rows;
		star.cellWidth = star.width/columns;
		star.cellHeight = star.height/rows;
		star.cellCount = columns*rows;

		while(star.cells.length < star.cellCount) {
			star.cells.push(new Cell());
		}

		for(y=rows-1; y>=0; y--) {
			for(x=columns-1; x>=0; x--) {
				star.cells[y*columns+x].reset(
					x*star.cellWidth,
					y*star.cellHeight,
					star.cellWidth,
					star.cellHeight
				);
				// cell = new Cell(
				// 	x*star.cellWidth,
				// 	y*star.cellHeight,
				// 	star.cellWidth,
				// 	star.cellHeight
				// );
				// star.cells.push(cell);
			}
		}

		// for(y=0; y<rows; y++) {
		// 	for(x=0; x<columns; x++) {
		// 		cell = star.cells[y*columns + x];
		// 		if(x > 0) {
		// 			cell.westernCell = star.cells[y*columns + x-1];
		// 		}
		// 		if(x < columns-1) {
		// 			cell.easternCell = star.cells[y*columns + x+1];
		// 		}
		// 		if(y > 0) {
		// 			cell.northernCell = star.cells[(y-1)*columns + x];
		// 		}
		// 		if(y < rows-1) {
		// 			cell.southernCell = star.cells[(y+1)*columns + x];
		// 		}
		// 	}
		// }

		star.ui.setGridSize(columns, rows);

		star.gridCellSize = star.cellWidth*star.cellHeight;

	},

	//
	//  Grid Resizing
	//_________________//

	previousGridCellSizes: [],

	previousSimTimes: [],

	defaultCellSize: 10000,

	gridSizingInterval: null,

	GRID_RESIZE_INTERVAL: 50,
	//#replace star.GRID_RESIZE_INTERVAL 50

	GRID_SIZING_VARIANCE: 0.1,
	//#replace star.GRID_SIZING_VARIANCE 0.1

	PREVIOUS_GRID_CELL_SIZES_MAX_SIZE: 20,
	//#replace star.PREVIOUS_GRID_CELL_SIZES_MAX_SIZE 20

	setAutomaticGridSizing: function star_setAutomaticGridSizing(on) {
		if(on && !star.enableAutomaticGridSizing) {
			star.startGridSizingRoutine();
		}
		else if(!on && star.enableAutomaticGridSizing) {
			star.clearGridSizingRoutine();
		}
		star.enableAutomaticGridSizing = on;
	},

	startGridSizingRoutine: function star_startGridSizingRoutine() {

		star.gridSizingInterval = setInterval(function() {

			var min = Infinity,
				minIndex,
				i, time,
				cellSize = star.defaultCellSize;

			if(star.previousGridCellSizes.length) {
				for(i=star.previousGridCellSizes.length-1; i>=0; i--) {
					time = star.previousSimTimes[i];
					if(time < min) {
						cellSize = star.previousGridCellSizes[i];
						min = time;
					}
				}
			}

			if(Math.random() < 0.5) {
				cellSize += cellSize * star.GRID_SIZING_VARIANCE;
			}
			else {
				cellSize -= cellSize * star.GRID_SIZING_VARIANCE;
			}

			star.setGridCellSize(cellSize);

		}, star.GRID_RESIZE_INTERVAL);

	},

	clearGridSizingRoutine: function star_clearGridSizingRoutine() {
		if(star.gridSizingInterval !== null) {
			clearInterval(star.gridSizingInterval);
		}
	},

	setGridCellSize: function star_setGridCellSize(size) {

		var canvasSize = star.width * star.height,
			count = canvasSize / size,
			actualSize = canvasSize / count;

		// var columns = Math.round(star.width / Math.sqrt(actualSize));
		// var rows = Math.round(star.height / Math.sqrt(actualSize));

		star.setGridSize(
			Math.round(star.width / Math.sqrt(actualSize)),
			Math.round(star.height / Math.sqrt(actualSize))
		);

	},

	//
	//  Physics
	//___________//

	moveForceAll: function star_moveForceAll(dt) {

		if(isNaN(dt)) {
			dt = 0.01;
		}

		dt *= star.gravity;

		var bodies = star.bodies,
			cells = star.cells,
			i, body,
			cell, cellBodies,
			i2,
			xDiff, yDiff,
			qdist, force,
			cellMass, cellMassX, cellMassY;

		// Set the force on each body
		// for( i=bodies.length-1; i>=0; i--) {
		// 	body = bodies[i];
		// 	// Drop forces
		// 	body.fx = 0;
		// 	body.fy = 0;
		// }

		// Loop through each cell
		for( i=star.cellCount-1; i>=0; i--) {
			cell = cells[i];
			cellBodies = cell.bodies;
			cellMass = cell.mass;

			// Loop over each body in the cell
			for( i2=cellBodies.length-1; i2>=0; i2-- ) {

				// Apply the gravity of every other body in the cell
				star.gravi(cellBodies[i2], cellBodies);

				// And nearby cells
				// if(cell.southernCell !== null) {
				// 	star.gravi(cellBodies[i2], cell.southernCell.bodies);
				// }
				// if(cell.northernCell !== null) {
				// 	star.gravi(cellBodies[i2], cell.northernCell.bodies);
				// }
				// if(cell.westernCell !== null) {
				// 	star.gravi(cellBodies[i2], cell.westernCell.bodies);
				// }
				// if(cell.easternCell !== null) {
				// 	star.gravi(cellBodies[i2], cell.easternCell.bodies);
				// }

			}

			if(cellMass) {
				cellMassX = cell.massX;
				cellMassY = cell.massY;
				for( i2=bodies.length-1; i2>=0; i2--) {
					body = bodies[i2];
					if(body.cell !== cell) {

						xDiff = body.x - cellMassX;
						yDiff = body.y - cellMassY;

						qdist = xDiff*xDiff + yDiff*yDiff;

						// qdist = (body.mass * cellMass / qdist)/Math.sqrt(qdist);

						// body.fx -= xDiff*qdist;
						// body.fy -= yDiff*qdist;

						force = 2 * body.mass * cellMass / qdist;

						body.fx -= force*(body.x - cell.massX)/Math.sqrt(qdist);
						body.fy -= force*(body.y - cell.massY)/Math.sqrt(qdist);

					}
				}
			}

		}

		// Apply the force on each body
		for( i=bodies.length-1; i>=0; i--) {
			body = bodies[i];
			if(body.canMove) {
				body.applyForce(dt);
			}
		}

	},

	gravi: function star_gravi(body, bodies) {

		var bx = body.x,
			by = body.y,
			bmass = body.mass,
			xDiff, yDiff,
			i, current_body, qdist, force;

		for( i=bodies.length-1; i>=0; i--) {
			current_body = bodies[i];
			if(body !== current_body) {

				xDiff = bx - current_body.x;
				yDiff = by - current_body.y;

				// Only works if not at equivalent points
				if(xDiff && yDiff) {

					qdist = xDiff*xDiff + yDiff*yDiff;

					// force = ( bmass * current_body.mass / qdist ) / Math.sqrt(qdist);

					// body.fx -= xDiff*force;
					// body.fy -= yDiff*force;

					force = 2 * bmass * current_body.mass / qdist;

					body.fx -= force*(body.x - current_body.x)/Math.sqrt(qdist);
					body.fy -= force*(body.y - current_body.y)/Math.sqrt(qdist);

				}

			}
		}

	},

	// cellToCellCollisions: function star_cellToCellCollisions(bodies1, bodies2) {

	// 	return;

	// 	var i, j,
	// 		b1, b1x, b1y, b1r,
	// 		b2, b2r,
	// 		rx, ry,
	// 		alpha, c, proj, proj_rx, proj_ry, delta,
	// 		xDiff, yDiff, dist_sq,
	// 		energy, index,
	// 		bounce = star.particleBounce,
	// 		temp;

	// 	if(star.flipNextCollisionCheckDirection) {
	// 		temp = bodies1;
	// 		bodies2 = bodies1;
	// 		bodies1 = temp;
	// 	}

	// 	// Loop over every body in cell 1
	// 	// for( i=0; i<bodies1.length; i++ ) {
	// 	for( i=bodies1.length-1; i>=0; i-- ) {

	// 		b1 = bodies1[i];

	// 		if(b1.canMove) {

	// 			b1x = b1.x;
	// 			b1y = b1.y;
	// 			b1r = b1.radius;

	// 			// Loop over every body in cell 2
	// 			// for(j=0; j<bodies2.length; j++) {
	// 			for( j=bodies2.length-1; j>=0; j-- ) {

	// 				b2 = bodies2[j];

	// 				if(b2.canMove) {

	// 					b2r = b2.radius;

	// 					xDiff = b1x-b2.x;
	// 					yDiff = b1y-b2.y;
	// 					dist_sq = xDiff*xDiff + yDiff*yDiff;

	// 					// Ensure that b2 can move and that the distance is not 0, yet within the radius
	// 					if(dist_sq !== 0 && dist_sq < (b1r+b2r)*(b1r+b2r)) {

	// 						energy = b1.mass*(b1.vx*b1.vx + b1.vy*b1.vy) + b2.mass*(b2.vx*b2.vx + b2.vy*b2.vy);

	// 						if(energy > b1.threshold+b2.threshold && b1.type+b2.type <= 7) {
	// 							index = star.bodies.indexOf(b2);
	// 							star.bodies.splice(index, 1);
	// 							b1.setType(b1.type+b2.type);
	// 							b1.vx /= 4;
	// 							b1.vy /= 4;
	// 							continue;
	// 						}

	// 						alpha = (b1.mass/b2.mass)/2 + 0.5;
	// 						c = Math.sqrt(dist_sq); //normalize vectors

	// 						//normalized radius - vector
	// 						rx = (b2.x-b1x)/c;
	// 						ry = (b2.y-b1y)/c;

	// 						proj = (b2.vx-b1.vx) * rx + (b2.vy-b1.vy) * ry;

	// 						proj_rx = proj * rx;
	// 						proj_ry = proj * ry;
							
	// 						b2.vx += (bounce*proj_rx * (1 - alpha) / alpha) - proj_rx;
	// 						b2.vy += (bounce*proj_ry * (1 - alpha) / alpha) - proj_ry;

	// 						b1.vx += (bounce*proj_rx) / alpha;
	// 						b1.vy += (bounce*proj_ry) / alpha;
							
	// 						delta = (b1r + b2r - c)/2;

	// 						b2.x += delta * rx;
	// 						b2.y += delta * ry;
							
	// 						b1x -= delta * rx;
	// 						b1y -= delta * ry;

	// 					}
	// 				}
	// 			}

	// 			b1.x = b1x;
	// 			b1.y = b1y;

	// 		}
	// 	}

	// },

	collisions: function star_collisions(bodies) {

		var i, j,
			b1, b1x, b1y, b1r,
			b2,
			rx, ry,
			alpha, c, proj, proj_rx, proj_ry, delta,
			xDiff, yDiff, dist_sq,
			energy, index,
			bounce = star.particleBounce;

		for( i=bodies.length-1; i>=0; i-- ) {

			b1 = bodies[i];

			if(b1.canMove) {

				b1x = b1.x;
				b1y = b1.y;
				b1r = b1.radius;

				for(j=i+1; j<bodies.length; j++) {

					b2 = bodies[j];

					if(b2.canMove) {

						xDiff = b1x-b2.x;
						yDiff = b1y-b2.y;
						dist_sq = xDiff*xDiff + yDiff*yDiff;

						if(dist_sq !== 0 && dist_sq < (b1r+b2.radius)*(b1r+b2.radius)) {

							energy = b1.mass*(b1.vx*b1.vx + b1.vy*b1.vy) + b2.mass*(b2.vx*b2.vx + b2.vy*b2.vy);

							if(energy > b1.threshold+b2.threshold && b1.type+b2.type <= 7) {
								index = star.bodies.indexOf(b2);
								star.bodies.splice(index, 1);
								star.classes[b1.type].count--;
								star.classes[b2.type].count--;
								star.classes[b1.type+b2.type].count++;
								b1.setType(b1.type+b2.type);
								star.ui.updateClassCounts();
								b1.vx *= b2.mass / b1.mass;
								b1.vy *= b2.mass / b1.mass;
								continue;
							}

							alpha = (b1.mass/b2.mass)/2 + 0.5;
							c = Math.sqrt(dist_sq); //normalize vectors

							//normalized radius - vector
							rx = (b2.x-b1x)/c;
							ry = (b2.y-b1y)/c;

							proj = (b2.vx-b1.vx) * rx + (b2.vy-b1.vy) * ry;

							proj_rx = proj * rx;
							proj_ry = proj * ry;
							
							b2.vx += (bounce*proj_rx * (1 - alpha) / alpha) - proj_rx;
							b2.vy += (bounce*proj_ry * (1 - alpha) / alpha) - proj_ry;

							b1.vx += (bounce*proj_rx) / alpha;
							b1.vy += (bounce*proj_ry) / alpha;
							
							delta = (b1r + b2.radius - c)/2;

							b2.x += delta * rx;
							b2.y += delta * ry;

							b1x -= delta * rx;
							b1y -= delta * ry;

						}
					}
				}

				b1.x = b1x;
				b1.y = b1y;

			}
		}

	},

	moveConstrainsAll: function star_moveConstrainsAll() {

		var cells = star.cells,
			i, cell,
			x, y,
			rows = star.rowCount,
			columns = star.columnCount;

		// Constrain the bodies to the field and sort them into the grid
		star.constrainAndSortBodies();

		// Handle collisions of bodies within each cell
		for( i=star.cellCount-1; i>=0; i-- ) {
			star.collisions(cells[i].bodies);
		}

		// Next we check for collisions between adjacent cells, overlapping collision checks between adjacent cells
		// This takes care of the possibility that two particles right on the edge of a cell are assured to collide

		// for( y=0; y<rows; y++ ) {
		// for( y=rows-1; y>=0; y-- ) {

		// 	// Overlap across the x dimension
		// 	// for( x=0; x<columns-1; x++ ) {
		// 	for( x=columns-2; x>=0; x-- ) {
		// 		star.cellToCellCollisions(cells[y*columns+x].bodies, cells[y*columns+x+1].bodies);
		// 	}

		// 	if(y < rows-1) {
		// 		// Overlap across the y dimension
		// 		// for( x=0; x<columns; x++ ) {
		// 		for( x=columns-1; x>=0; x-- ) {
		// 			star.cellToCellCollisions(cells[y*columns+x].bodies, cells[(y+1)*columns+x].bodies);
		// 		}
		// 	}

		// }

		// FIXME: Well dang, flipping the direction between frames causes the shakky particle problem
		// star.flipNextCollisionCheckDirection = !star.flipNextCollisionCheckDirection;

	},

	render: function star_render() {

		var bodies = star.bodies,
			cells,
			ctx = star.ctx,
			i, cell, body,
			lastClass = NaN;

		// Clear the canvas
		ctx.clearRect(0, 0, star.width, star.height);

		// Draw the grid if enabled
		if(star.drawGrid) {
			cells = star.cells;
			ctx.strokeStyle = 'white';
			ctx.fillStyle = 'red';
			for( i=star.cellCount-1; i>=0; i--) {
				cell = cells[i];
				ctx.strokeRect(cell.x, cell.y, star.cellWidth, star.cellHeight);
				if(star.drawCellMass) {
					ctx.fillText(cell.mass+'', cell.x+2, cell.y+10);
				}
				if(star.drawCellCenterOfMass && cell.mass) {
					ctx.fillRect(cell.massX-2, cell.massY-2, 4, 4);
				}
			}
		}

		// Loop through every body
		for( i=bodies.length-1; i>=0; i--) {
			body = bodies[i];
			// Set the color, unless its the same as the previously applied color
			if(lastClass !== body.type) {
				ctx.fillStyle = star.classes[body.type].color;
				lastClass = body.type;
			}
			// Draw the particle as a full arc
			ctx.beginPath();
			ctx.arc(body.x, body.y, body.radius * star.particleDrawSizeMultiplier, 0, 2*Math.PI, true);
			ctx.fill();
		}

	},

	updateParticles: function star_updateParticles() {
		for( var i=star.bodies.length-1; i>=0; i--) {
			star.bodies[i].setType(star.bodies[i].type);
		}
	},

	//
	//  User Interface
	//__________________//

	ui: {

		init: function star_ui_init() {

			this.cpMain = document.getElementById('cp-main');

			var headerButtons = document.getElementsByClassName('header-button');
			for( var i=0; i<headerButtons.length; i++ ) {
				headerButtons[i].onclick = function(e) {
					var target = $('#'+e.target.getAttribute('data-toggle'));
					if(target.is(':visible')) {
						target.hide();
					}
					else {
						target.show();
					}
				};
			}

			// Main

			var restart = $('#restart');
			restart.click(star.restart);

			var pause = $('#pause');
			pause.click(function() {
				if(star.running) {
					star.stop();
					pause.text('Resume');
				}
				else {
					star.start();
					pause.text('Pause');
				}
			});

			$('#cp-main-count').val(star.particleCount);
			$('#cp-main-count').change(function(e) {
				star.particleCount = parseInt($(e.target).val(), 10);
			});

			// Physics

			$('#cp-physics-gravity').val(star.gravity);
			$('#cp-physics-gravity').change(function(e) {
				star.gravity = parseFloat($(e.target).val());
			});

			$('#cp-physics-particleBounce').val(star.particleBounce);
			$('#cp-physics-particleBounce').change(function(e) {
				star.particleBounce = parseFloat($(e.target).val());
			});

			$('#cp-physics-wallBounce').val(star.wallBounce);
			$('#cp-physics-wallBounce').change(function(e) {
				star.wallBounce = parseFloat($(e.target).val());
			});

			$('#cp-physics-wallFrict').val(star.wallFrict);
			$('#cp-physics-wallFrict').change(function(e) {
				star.wallFrict = parseFloat($(e.target).val());
			});

			// Grid

			var enableAutomaticGridSizing = $('#enableAutomaticGridSizing');
			enableAutomaticGridSizing.attr('checked', star.enableAutomaticGridSizing);
			if(star.enableAutomaticGridSizing) {
				$('#cp-gridSize-columns').attr('disabled', true);
				$('#cp-gridSize-rows').attr('disabled', true);
			}
			enableAutomaticGridSizing.change(function(e) {
				var on = e.target.checked;
				star.setAutomaticGridSizing(on);
				if(on) {
					$('#cp-gridSize-columns').attr('disabled', true);
					$('#cp-gridSize-rows').attr('disabled', true);
					star.ui.setGridSize(star.columnCount, star.rowCount);
				}
				else {
					$('#cp-gridSize-columns').attr('disabled', false);
					$('#cp-gridSize-rows').attr('disabled', false);
				}
			});

			$('#cp-gridSize-columns').change(function(e) {
				star.setGridSize(parseInt($(e.target).val(), 10), star.rowCount);
			});

			$('#cp-gridSize-rows').change(function(e) {
				star.setGridSize(star.columnCount, parseInt($(e.target).val(), 10));
			});

			// Rendering

			$('#cp-drawGrid')[0].checked = star.drawGrid;
			$('#cp-drawGrid').change(function(e) {
				star.drawGrid = e.target.checked;
			});

			$('#cp-drawCellMass')[0].checked = star.drawCellMass;
			$('#cp-drawCellMass').change(function(e) {
				star.drawCellMass = e.target.checked;
			});

			$('#cp-drawCellCenterOfMass')[0].checked = star.drawCellCenterOfMass;
			$('#cp-drawCellCenterOfMass').change(function(e) {
				star.drawCellCenterOfMass = e.target.checked;
			});

			$('#cp-rendering-particleDrawSizeMultiplier').val(star.particleDrawSizeMultiplier);
			$('#cp-rendering-particleDrawSizeMultiplier').change(function(e) {
				star.particleDrawSizeMultiplier = parseFloat($(e.target).val());
			});

			setInterval(function() {
				$('#simTime').text(star.simTime);
				$('#drawTime').text(star.drawTime);
			}, 500);

			// Particle Classes
			var data, elem;
			for(i=0;i<star.classes.length-1;i++) {

				data = star.classes[i];

				elem = $('#cp-class'+i);
				if(elem.length) {

					elem = $('#cp-class'+i+'-color');
					elem[0].data = data;
					elem.val(data.color);
					elem.change(function(e) {
						e.target.data.color = $(e.target).val();
						// star.updateParticles();
					});

					elem = $('#cp-class'+i+'-physicalDiameter');
					elem[0].data = data;
					elem.val(data.diameter);
					elem.change(function(e) {
						e.target.data.diameter = parseInt($(e.target).val(), 10);
						star.updateParticles();
					});

					elem = $('#cp-class'+i+'-mass');
					elem[0].data = data;
					elem.val(data.mass);
					elem.change(function(e) {
						e.target.data.mass = parseInt($(e.target).val(), 10);
						star.updateParticles();
					});

					elem = $('#cp-class'+i+'-threshold');
					elem[0].data = data;
					elem.val(data.threshold);
					elem.change(function(e) {
						e.target.data.threshold = parseInt($(e.target).val(), 10);
						star.updateParticles();
					});


				}

			}

		},

		setGridSize: function star_ui_setGridSize(columns, rows) {
			$('#cp-gridSize-columns').val(columns);
			$('#cp-gridSize-rows').val(rows);
		},

		updateClassCounts: function star_ui_updateClassCounts() {
			for(var i=1;i<star.classes.length-1;i++) {
				$('#class'+i+'_count').text(star.classes[i].count);
			}
		}

	}

};

})();