.PHONY: all

all:
	bash raytrace_all.sh
	convert -delay 50 -loop 10 *.png animation.gif

update:
	cp ~/PhD/projects/pov-ray/library.inc .
	cp ~/PhD/projects/pov-ray/scene.inc .

clean:
	rm *.pov
	rm *.gif
	rm *.png
