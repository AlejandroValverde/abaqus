#Variables, clean

# - N: Number of half-spanwise elements
# - width_element: Width of the elements in the spanwise direction
# - y, y_C: Control points in the spanwise direction
# - c_0, c: Chord length [m]
# - span: Wing span [m]
# - width: width of the wing (half of the span width) [m]
# - PositionFrontSpar: Dimensional-less (c) position of the front spar
# - PositionRearSpar: Dimensional-less (c) position of the rear spar
# - x_spar1: PositionFrontSpar*c
# - x_spar2: PositionRearSpar*c
# - anzRibs: 
# - s_stringer: Stringer thickness
# - element_width: width/N
# - x_str1, x_str2, x_str3: Position of str on chord 1, 2, 3
# - x_spar2: 
# - n_skinpts: Number of points in the skin
# - skin_points_upX: Point vector for the skin 
# - skin_points_upY: Point vector for the skin
# - skin_points_dnX: Point vector for the skin
# - skin_points_dnY: Point vector for the skin
# - mod6_UP_r, mod6_UP_m, mod6_UP_f: Define nodes to find edges in the modeling part of the code, upper surface
# - mod6_DN_m: Define nodes to find edges in the modeling part of the code, lower surface
# - mod7_UP, mod7_DN:
# - mod8_UP, mod8_DN: 
# - x_surf, y_surf: Surface points
# - points_nose, points_box_up, points_box_low, points_trail1_low, points_trail1_up, 
#	points_trail2_low, points_trail2_up, points_trail3_low, points_trail3_up, 
#	points_trail4_low, points_trail4_up: Points defining sections of the airfoil p.e.: x=points_nose[i][0], ypoints_nose[i][1]
# - points_trail_rib_low, points_trail_rib_up: Points defining sections of the rib
# - rib_box_up, rib_box_low: 
# 