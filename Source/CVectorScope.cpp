
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"

namespace Signalizer
{
	// swapping the right channel might give an more intuitive view
	#define SWAP_RIGHT

	enum RenderButtons
	{
		fillLines,
		antiAlias,
		fadeHistory,
		rgbHistory
	};

	std::unique_ptr<juce::Component> CVectorScope::createEditor()
	{
		auto content = new Signalizer::CContentPage();

		if (auto page = content->addPage("Settings", "icons/svg/wrench.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&ktransform, 0);
				page->addSection(section, "Transform");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kgain, 0);
				section->addControl(&kwindow, 1);
				page->addSection(section, "Utility");
			}
		}

		if (auto page = content->addPage("Rendering", "icons/svg/painting.svg"))
		{
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kantiAlias, 0);
				section->addControl(&kfadeOld, 1);
				section->addControl(&kdrawLines, 2);
				page->addSection(section, "Options");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kdrawingColour, 0);
				section->addControl(&kgraphColour, 0);
				section->addControl(&kbackgroundColour, 0);
				section->addControl(&kskeletonColour, 0);
				section->addControl(&krotation, 1);
				section->addControl(&kzrotation, 1);
				section->addControl(&kdrawGraph, 1);
				page->addSection(section, "Look");
			}
		}
		return std::unique_ptr<juce::Component>(content);

	}

	CVectorScope::CVectorScope(cpl::AudioBuffer & data)
	:
		audioStream(data),
		windowSize(1000),
		kwindow("Window Size", cpl::CKnobSlider::ControlType::ms),
		krotation("Matrix Rotation"),
		kgain("Gain"),
		kgraphColour("Graph colour"),
		kbackgroundColour("Background colour"),
		kdrawingColour("Drawing colour"),
		kskeletonColour("Skeleton colour"),
		kzrotation("z-plane"),
		processorSpeed(0), 
		audioStreamCopy(2),
		lastFrameTick(0), 
		isFrozen(false),
		lastMousePos()
	{
		setOpaque(true);
		textbuf = std::unique_ptr<char>(new char[300]);
		processorSpeed = juce::SystemStats::getCpuSpeedInMegaherz();
		initPanelAndControls();
	}

	void CVectorScope::suspend()
	{


	}

	void CVectorScope::resume()
	{


	}

	juce::Component * CVectorScope::getWindow()
	{
		return this;
	}

	CVectorScope::~CVectorScope()
	{
		notifyDestruction();
	}

	void CVectorScope::initPanelAndControls()
	{
		// listeners
		kwindow.bAddPassiveChangeListener(this);
		kwindow.bAddFormatter(this);
		kgain.bAddFormatter(this);
		krotation.bAddPassiveChangeListener(this);
		krotation.bAddFormatter(this);
		kzrotation.bAddFormatter(this);
		ktransform.bAddPassiveChangeListener(this);
		// buttons
		kantiAlias.bSetTitle("Antialias lines");
		kantiAlias.setToggleable(true);
		kfadeOld.bSetTitle("Fade older points");
		kfadeOld.setToggleable(true);
		kdrawLines.bSetTitle("Draw lines");
		kdrawLines.setToggleable(true);
		kdrawGraph.bSetTitle("Draw graph");
		kdrawGraph.setToggleable(true);
		// descriptions.
		kwindow.bSetDescription("The size of the displayed time window.");
		kgain.bSetDescription("How much the input is scaled (or the input gain), effectively zooms the view.");
		krotation.bSetDescription("The amount of degrees to rotate the x/y points around the origin.");
		kantiAlias.bSetDescription("Antialiases lines (if set). May slow down rendering.");
		kfadeOld.bSetDescription("If set, gradually older samples will be faded linearly.");
		kdrawLines.bSetDescription("If set, interconnect samples linearly.");
		kdrawingColour.bSetDescription("The main colour to paint with.");
		kgraphColour.bSetDescription("The colour of the graph.");
		kbackgroundColour.bSetDescription("The background colour of the view");
		kdrawGraph.bSetDescription("Draws a graph in the view, if enabled.");
		kskeletonColour.bSetDescription("The colour of the box skeleton indicating the OpenGL camera clip box.");
		// design
		setMouseCursor(juce::MouseCursor::DraggingHandCursor);
	}

	void CVectorScope::save(cpl::CSerializer::Archiver & archive, long long int version)
	{
		archive << kwindow;
		archive << kgain;
		archive << krotation;
		archive << kantiAlias;
		archive << kfadeOld;
		archive << kdrawGraph;
		archive << kdrawLines;
		archive << kgraphColour;
		archive << kbackgroundColour;
		archive << kdrawingColour;
		archive << ktransform.getTransform3D();
		archive << kskeletonColour;
	}

	void CVectorScope::load(cpl::CSerializer::Builder & builder, long long int version)
	{
		builder >> kwindow;
		builder >> kgain;
		builder >> krotation;
		builder >> kantiAlias;
		builder >> kfadeOld;
		builder >> kdrawGraph;
		builder >> kdrawLines;
		builder >> kgraphColour;
		builder >> kbackgroundColour;
		builder >> kdrawingColour;
		builder >> ktransform.getTransform3D();
		builder >> kskeletonColour;


		ktransform.syncEditor();

	}

	void CVectorScope::paint(juce::Graphics & g)
	{
		auto cStart = cpl::Misc::ClockCounter();
		auto colour = juce::Colour(kgraphColour.getControlColourAsColour());
		
		auto & matrix = ktransform.getTransform3D();
		
		auto currentXOffset = (matrix.position.x) * getWidth() / 2.0f;
		auto currentYOffset = (matrix.position.y) * -getHeight() / 2.0f; // invert

		if (kdrawGraph.bGetValue() > 0.5)
		{
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());
			auto middleHeight = getHeight() / 2.0f;
			auto middleWidth = getWidth() / 2.0f;


			g.setColour(colour);

			g.drawLine(0.f, middleHeight + currentYOffset, (float)getWidth(), middleHeight + currentYOffset);
			g.drawLine(middleWidth + currentXOffset, 0.f, middleWidth + currentXOffset, (float)getHeight());

			g.addTransform(juce::AffineTransform::translation(currentXOffset, currentYOffset));



			#ifdef SWAP_RIGHT
				static const char * channels[] = { "+L", "+R", "-L", "-R" };
			#else
				static const char * channels[] = { "+L", "-R", "-L", "+R" };
			#endif

			const auto rotation = krotation.bGetValue() * 2 * M_PI;
			for (int i = 0; i < 4; ++i)
			{
				auto ycoord = sin(i * M_PI * 0.5 + -rotation);
				auto xcoord = sin(i * M_PI * 0.5 + -rotation + M_PI / 2);

				g.drawSingleLineText(channels[i],
					cpl::Math::round<int>(middleWidth * 0.05 + (xcoord * middleWidth + middleWidth) * 0.9), // make the ellipsis a bit smaller to make space for text
					cpl::Math::round<int>(middleHeight * 0.05 + (ycoord * middleHeight + middleHeight) * 0.9));

			}

			// revert
			g.addTransform(juce::AffineTransform::translation(-currentXOffset, -currentYOffset));
			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu", 
				getWidth(), getHeight(), fps, cpuTime);
			g.drawSingleLineText(textbuf.get(), 10, 20);
			
		}
	}

	void CVectorScope::renderOpenGL()
	{
		if (audioStream.empty())
			return;

		auto cStart = cpl::Misc::ClockCounter();
		OpenGLHelpers::clear(kbackgroundColour.getControlColourAsColour());

		const bool fillPath = kdrawLines.bGetValue() > 0.5;
		const bool fadeHistory = kfadeOld.bGetValue() > 0.5;
		const bool antiAlias = kantiAlias.bGetValue() > 0.5;

		const std::size_t numSamples = audioStream[0].size - 1;

		cpl::AudioBuffer * buffer;
		juce::Colour colour;
		float red = 0, blue = 0, green = 0, alpha = 0;
		auto loadColour = [&](juce::Colour colour)
		{
			red = colour.getFloatRed(); green = colour.getFloatGreen(); blue = colour.getFloatBlue(); alpha = colour.getFloatAlpha();
			glColor4f(red, green, blue, alpha);
		};

		if (isFrozen)
		{

			buffer = &audioStreamCopy;
		}
		else if (isSynced)
		{
			std::vector<cpl::CMutex> locks;

			for (auto & buffer : audioStream)
				locks.emplace_back(buffer);

			for (auto i = 0; i < audioStream.size(); ++i)
			{
				audioStream[i].clone(audioStreamCopy[i]);

			}
			buffer = &audioStreamCopy;
		}
		else
		{
			buffer = &audioStream;

		}

		/*	http://www.juce.com/forum/topic/how-do-i-draw-text-over-openglcomponent?page=1
			Set up rendering
		*/
		float xn, yn, sleft, sright;
		float gain = getGain();
		xn = yn = sleft = sright = 0.0f;
		float sampleFade = 1.0 / numSamples;
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_COLOR);

		loadColour(kdrawingColour.getControlColourAsColour());

		
		//glTranslatef(xoffset + xdrag, yoffset + ydrag, numSamples / audioStream[0].sampleRate);

		glLoadIdentity();
		glMatrixMode(GL_PROJECTION);
		auto & matrix = ktransform.getTransform3D();
		matrix.applyToOpenGL();
		glFrustum(-1, 1, -1, 1, 0, 5);
		//antiAlias ? glShadeModel(GL_SMOOTH) : glShadeModel(GL_FLAT);
		//antiAlias ? glEnable(GL_LINE_SMOOTH) : glDisable((GL_LINE_SMOOTH));
		//antiAlias ? glEnable(GL_MULTISAMPLE) : glDisable((GL_MULTISAMPLE));
		fillPath ? glBegin(GL_LINE_STRIP) : glBegin(GL_POINTS);

		cpl::CAudioBuffer::CChannelBuffer * leftChannel, *rightChannel;
		leftChannel = &((*buffer)[0]);
		rightChannel = &((*buffer)[1]);
		float fade = 0;
		float z = kzrotation.bGetValue();
		
		if (fadeHistory)
		{
			for (std::size_t i = 0; i < numSamples; ++i)
			{
				sleft = leftChannel->singleCheckAccess(i);
				#ifdef SWAP_RIGHT
					sright = -rightChannel->singleCheckAccess(i);
				#else
					sright = rightChannel->singleCheckAccess(i);
				#endif
				//fade = i * sampleFade * 0.5f;
				//glColor3f(fade * red, fade * green, fade * blue);

				xn = sleft * cosrol - sright * sinrol;
				yn = sleft * sinrol + sright * cosrol;
				glVertex3f(xn * gain, yn * gain, i * sampleFade - 1);
				//glVertex2d(xn * gain, yn * gain);

			}
			//if (numSamples & 0x1) // we need to supply an even pair of vertices when drawing lines
			//	glVertex2d(xn, yn);
		}
		else
		{
			for (std::size_t i = 0; i < numSamples; ++i)
			{
				// use glDrawArrays instead here
				sleft = leftChannel->singleCheckAccess(i);
				#ifdef SWAP_RIGHT
					sright = -rightChannel->singleCheckAccess(i);
				#else
					sright = rightChannel->singleCheckAccess(i);
				#endif
				// todo: check here for 0, 0 coordinates
				xn = sleft * cosrol - sright * sinrol;
				yn = sleft * sinrol + sright * cosrol;

				glVertex2f(xn * gain, yn * gain);

			}
		}
		glEnd();
		glBegin(GL_LINE_STRIP);

		glColor3f(1.f, 1.f, 1.f);
		
		glVertex3f(0.0f, 0.0f, 0.0f);  glVertex3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f); glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);	glVertex3f(0.0f, -1.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f); glVertex3f(0.0f, 0.0f, 0.0f);
		
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f); glVertex3f(0.0f, 0.0f, -1.0f);
		glVertex3f(0.0f, 0.0f, 0.0f); //glVertex3f(0.0f, 0.0f, 0.0f);
		/*
			Stop rendering
		*/
		glEnd();
		
		loadColour(kskeletonColour.getControlColourAsColour());
		glBegin(GL_LINES);
		// draw skeleton
		int nlines = 14;
		auto rel = 1.0f / nlines;

		// front
		for (int i = 0; i <= nlines; ++i)
		{
			glVertex3f(i * rel * 2 - 1, -1.f, 0.0f);
			glVertex3f(i * rel * 2 - 1, 1.f, 0.0f);
			//glVertex3f((i + 1) * rel * 2 - 1, 1.f, 0.0f);
		}

		for (int i = 0; i <= nlines; ++i)
		{
			glVertex3f(-1.0f, i * rel * 2 - 1, 0.0f);
			glVertex3f(1.0f, i * rel * 2 - 1, 0.0f);
			//glVertex3f(1.0f, (i + 1) * rel * 2 - 1, 0.0f);
		}
		// back
		for (int i = 0; i <= nlines; ++i)
		{
			glVertex3f(i * rel * 2 - 1, -1.f, -1.0f);
			glVertex3f(i * rel * 2 - 1, 1.f, -1.0f);
			//glVertex3f((i + 1) * rel * 2 - 1, 1.f, -1.0f);
		}

		for (int i = 0; i <= nlines; ++i)
		{
			glVertex3f(-1.0f, i * rel * 2 - 1, -1.0f);
			glVertex3f(1.0f, i * rel * 2 - 1, -1.0f);
			//glVertex3f(1.0f, (i + 1) * rel * 2 - 1, -1.0f);
		}

		glEnd();

		//glScalef(-0.25f, -0.25f, -0.25f);
		//glBegin(GL_LINE_STRIP);


		// face
		/*glVertex3f(0, 0, 0); glVertex3f(0, 0, 1); glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0); glVertex3f(1, 0, 1); glVertex3f(1, 0, 0);
		glVertex3f(1, 1, 0); glVertex3f(1, 1, 1); glVertex3f(1, 1, 0);
		glVertex3f(0, 1, 0); glVertex3f(0, 1, 1); glVertex3f(0, 1, 0);
		glVertex3f(0, 0, 0); glVertex3f(0, 0, 1); glVertex3f(0, 0, 0);

		glVertex3f(0, 0, 1);
		glVertex3f(1, 0, 1);
		glVertex3f(1, 1, 1);
		glVertex3f(0, 1, 1);
		glVertex3f(0, 0, 1);*/
		cpl::GraphicsND::Transform3D<float>::revert();

		renderCycles = cpl::Misc::ClockCounter() - cStart;

		auto tickNow = juce::Time::getHighResolutionTicks();
		avgFps.setNext(tickNow - lastFrameTick);
		lastFrameTick = tickNow;
	}

	void CVectorScope::freeze()
	{
		isFrozen = true;
		std::vector<cpl::CMutex> locks;
		// lock streams firstly if we are synced.
		if (isSynced)
		{
			for (auto & buffer : audioStream)
			{
				locks.emplace_back(buffer);
			}
		}

		for (unsigned i = 0; i < audioStream.size(); ++i)
		{
			audioStream[i].clone(audioStreamCopy[i]);
		}
	}

	void CVectorScope::unfreeze()
	{
		isFrozen = false;
	}


	void CVectorScope::resized()
	{


	}
	void CVectorScope::mouseWheelMove(const MouseEvent& event, const MouseWheelDetails& wheel)
	{
		auto amount = wheel.deltaY;
		auto normalizedSign = amount > 0 ? 1.0f : -1.0f;
		if (event.mods.isCtrlDown())
		{
			// increase gain
			kgain.bSetValue(kgain.bGetValue() + amount / 20);
		}
		else /* zoom graph */
		{
			
			auto mouseXCoord = (float(event.x) / getWidth()) * 2 - 1;
			auto mouseYCoord = (float(event.y) / getHeight()) * 2 - 1;

			auto & matrix = ktransform.getTransform3D();
			//matrix.scale.x += amount;
			//matrix.scale.y += amount;
			auto actualAmount = (1 + amount / 5) * matrix.scale.z;
			auto deltaIncrease = (actualAmount - matrix.scale.z) / matrix.scale.z;
			matrix.scale.z = actualAmount;

			matrix.scale.x *= 1 + deltaIncrease;
			matrix.scale.y *= 1 + deltaIncrease;
			//matrix.position.x *=  1 + normalizedSign * deltaIncrease;
			//matrix.position.y *= 1 + normalizedSign * deltaIncrease;
			ktransform.syncEditor();

		}
		
	}
	void CVectorScope::mouseDoubleClick(const MouseEvent& event)
	{
		
		if (event.mods.isLeftButtonDown())
		{
			// reset all zooming, offsets etc. when doubleclicking left
			kgain.bSetValue(0.5f);
			auto & matrix = ktransform.getTransform3D();
			matrix.position.y = matrix.position.x = 0;
			ktransform.syncEditor();
		}
	}
	void CVectorScope::mouseDrag(const MouseEvent& event)
	{
		auto & matrix = ktransform.getTransform3D();
		auto factor = float(getWidth()) / getHeight();
		auto deltaDifference = event.position - lastMousePos;
		if (event.mods.isCtrlDown())
		{
			
			auto deltaDifference = event.position - lastMousePos;
			
			matrix.rotation.y = std::fmod(deltaDifference.x * 0.3f + matrix.rotation.y, 360.f);
			
			matrix.rotation.x = std::fmod(deltaDifference.y * 0.3f + matrix.rotation.x, 360.f);
			//matrix.rotation.x = quat.vector.x * 180;
			//matrix.rotation.y = quat.vector.y * 180;
			//matrix.rotation.z = quat.vector.z * 90;


		}
		else
		{
			// needed to balance opengl coordinate system
			//ydrag = (-event.getDistanceFromDragStartY() / 500.0f) * factor;
			matrix.position.x += deltaDifference.x / 500.f;
			matrix.position.y += factor * -deltaDifference.y / 500.f;
		}
		ktransform.syncEditor();
		
		
		lastMousePos = event.position;
	}
	void CVectorScope::mouseUp(const MouseEvent& event)
	{
		if (event.mods.testFlags(ModifierKeys::rightButtonModifier))
		{
			unfreeze();
		}
	}
	void CVectorScope::mouseDown(const MouseEvent& event)
	{
		if (event.mods.testFlags(ModifierKeys::rightButtonModifier))
		{
			freeze();
		}
		lastMousePos = event.position;
	}

	bool CVectorScope::valueToString(const cpl::CBaseControl * ctrl, std::string & buffer, cpl::iCtrlPrec_t value)
	{
		char buf[200];
		if (ctrl == &kgain)
		{
			sprintf(buf, "%.2f dB (%d%%)", cpl::Math::fractionToDB(getGain()), (int)cpl::Misc::Round(getGain() * 100));
			buffer = buf;
			return true;
		}
		else if (ctrl == &kwindow)
		{
			auto bufLength = cpl::Math::round<int>( value * 1000);
			sprintf(buf, "%d ms", bufLength);
			buffer = buf;
			return true;
		}
		else if (ctrl == &krotation || ctrl == &kzrotation)
		{
			sprintf(buf, "%.2f degs", value * 360);
			buffer = buf;
			return true;

		}

		return false;
	}


	double CVectorScope::getGain()
	{
		return kgain.bGetValue() * 4;
		double val = 0;
		if (val > 0.5)
		{
			return 1 + (val * 2 - 1) * 9;

		}
		else if (val < 0.5)
		{
			return val * 2;
		}

		return 1;

	}

	double CVectorScope::mapScaleToFraction(double scale)
	{
		return scale / 4;
	}

	bool CVectorScope::stringToValue(const cpl::CBaseControl * ctrl, const std::string & buffer, cpl::iCtrlPrec_t & value)
	{

		if (ctrl == &kgain)
		{
			char * endPtr(nullptr);
			cpl::iCtrlPrec_t newVal = strtod(buffer.c_str(), &endPtr);
			if (endPtr > buffer.data())
			{
				value = mapScaleToFraction(cpl::Math::dbToFraction(newVal));
				return true;
			}
		}
		else if (ctrl == &kwindow)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(newValue / 1000.0, 0.0, 1.0);
				return true;
			}
		}
		else if (ctrl == &krotation || ctrl == &kzrotation)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(fmod(newValue, 360) / 360.0, 0.0, 1.0);
				return true;
			}
		}
		return false;
	}
	void CVectorScope::onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject)
	{
		// hmmm.....
		jassertfalse;
	}
	void CVectorScope::valueChanged(const cpl::CBaseControl * ctrl)
	{
		if (ctrl == &kwindow)
		{
			auto bufLength = ctrl->bGetValue() * 1000;
			for (auto & buffer : audioStream)
			{
				buffer.setLength(bufLength);
			}
			return;
		}
		else if (ctrl == &krotation)
		{
			double phase = ctrl->bGetValue() * 2 * M_PI;
			cosrol = (float)std::cos(phase);
			sinrol = (float)std::sin(phase);
			return;
		}	
		else if (ctrl == &ktransform)
		{
		}
	}


};
