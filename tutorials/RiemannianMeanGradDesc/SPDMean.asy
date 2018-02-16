	import three;
	import settings;

	surface ellipsoid(triple v1,triple v2,triple v3,real l1,real l2, real l3, triple pos=O) {
		transform3 T = identity(4);
		T[0][0] = l1*v1.x;
		T[1][0] = l1*v1.y;
		T[2][0] = l1*v1.z;
		T[0][1] = l2*v2.x;
		T[1][1] = l2*v2.y;
		T[2][1] = l2*v2.z;
		T[0][2] = l3*v3.x;
		T[1][2] = l3*v3.y;
		T[2][2] = l3*v3.z;
		T[0][3] = pos.x;
		T[1][3] = pos.y;
		T[2][3] = pos.z;
		return T*unitsphere;
	}

	size(200);

	real gDx=3.3584;
	real gDy=3.3584;
	real gDz=3.3584;

	currentprojection=perspective( camera=(gDx*0,gDy*0,1),up=Y, target=(gDx*0,gDy*0,gDz*0));
		draw(
			ellipsoid( (-0.20516,0.8418,0.49928), (0.35823,-0.41015,0.83872), (-0.91082,-0.35092,0.21741), 1.3683, 1.5432, 2.0782, (gDx*0,gDy*0,gDz*0)),rgb(0.59375,1,0));
