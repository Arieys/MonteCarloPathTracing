#include "morton code.h"

unsigned int expandBits(unsigned int v)
{
	v = (v * 0x00010001u) & 0xFF0000FFu;
	v = (v * 0x00000101u) & 0x0F00F00Fu;
	v = (v * 0x00000011u) & 0xC30C30C3u;
	v = (v * 0x00000005u) & 0x49249249u;
	return v;
}

unsigned int morton3D(float x, float y, float z)
{
	x = min(max(x * 1024.0f, 0.0f), 1023.0f);
	y = min(max(y * 1024.0f, 0.0f), 1023.0f);
	z = min(max(z * 1024.0f, 0.0f), 1023.0f);
	unsigned int xx = expandBits((unsigned int)x);
	unsigned int yy = expandBits((unsigned int)y);
	unsigned int zz = expandBits((unsigned int)z);
	return xx * 4 + yy * 2 + zz;
}

unsigned int getMortonCode(float x, float y, float z)
{
	float mx, my, mz;
	mx = (x - MINP) / (MAXP - MINP);
	my = (y - MINP) / (MAXP - MINP);
	mz = (z - MINP) / (MAXP - MINP);
	//map x,y,z to range [0,1]

	return morton3D(mx, my, mz);
}