import bpy
import math
import mathutils
import numpy as np
import sys, os
import json
#blender -b --python v1.py

#------------------------
# clear all
bpy.ops.wm.read_factory_settings(use_empty=True)
#------------------------

def make_path_curve(name, coords_list):
	# make a new curve
	crv = bpy.data.curves.new('crv', 'CURVE')
	crv.dimensions = '3D'

	# make a new spline in that curve
	spline = crv.splines.new(type='POLY')

	# a spline point for each point
	spline.points.add(len(coords_list)-1) # theres already one point by default

	# assign the point coordinates to the spline points
	for p, new_co in zip(spline.points, coords_list):
	    p.co = (new_co[0], new_co[1], new_co[2], 1) # (add nurbs weight)

	# make a new object with the curve
	obj = bpy.data.objects.new(name, crv)
	bpy.context.collection.objects.link(obj)

def cylinder_between(x1, y1, z1, x2, y2, z2, r):

	dx = x2 - x1
	dy = y2 - y1
	dz = z2 - z1    
	dist = math.sqrt(dx**2 + dy**2 + dz**2)

	bpy.ops.mesh.primitive_cylinder_add(
	  radius = r, 
	  depth = dist,
	  location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)   
	) 

	phi = math.atan2(dy, dx) 
	theta = math.acos(dz/dist) 

	bpy.context.object.rotation_euler[1] = theta 
	bpy.context.object.rotation_euler[2] = phi 

def make_space_time_rod(x, y, t):
	for i in range(0, len(t)-1):
		cylinder_between(x[i], y[i], t[i], x[i+1], y[i+1], t[i+1], float(0.25))


f = open("agents.json", "r")
scene = json.loads(f.read())
for rod in scene["agents"]:
	rod["v"] = np.array(rod["v"])

	x = rod["v"][:,0]
	y = rod["v"][:,1]
	t = rod["v"][:,2]


	h = np.linspace(t[0], t[-1], 10)
	hx = np.interp(h, t, x, left=None, right=None, period=None)
	hy = np.interp(h, t, y, left=None, right=None, period=None)
	hz = 0*hy;

	make_space_time_rod(x, y, t)
	make_path_curve("NurbsPath", np.array([hx, hy, hz], order='F').transpose())

	#add cube to act as pseudo car to animate along path
	imported_object = bpy.ops.import_scene.obj(filepath="scene_1/agent_1.obj", split_mode='OFF')
	obj_object = bpy.context.selected_objects[0] ####<--Fix
	bpy.context.view_layer.objects.active = obj_object

	#add constraint and specify properties
	bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN', center='MEDIAN')
	bpy.ops.object.constraint_add(type='FOLLOW_PATH')													#add follow path constraint
	bpy.data.objects[obj_object.name].constraints["Follow Path"].target = bpy.data.objects["NurbsPath"]			#sets the cube on the nurbs curve path
	bpy.data.objects[obj_object.name].constraints["Follow Path"].forward_axis = 'FORWARD_X'						#sets the face of the cube on the positive x-axis as front of car
	bpy.data.objects[obj_object.name].constraints["Follow Path"].use_curve_follow = True							#car will follow on path
	override={'constraint':bpy.data.objects[obj_object.name].constraints["Follow Path"]}
	bpy.ops.constraint.followpath_path_animate(override,constraint='Follow Path')

#animate the 'car'
#initial frame
# bpy.ops.object.select_all(action='DESELECT')														#move focus to path
# bpy.ops.object.select_name(name = "NurbsPath")

# bpy.data.scenes['Scene'].frame_current = 1															#set frame
# bpy.data.curves['NurbsPath'].eval_time = 0	

bpy.ops.wm.save_mainfile(filepath='./test.blend')

