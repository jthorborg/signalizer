
#include "CVectorScope.h"
#include <cstdint>
#include <cpl/CMutex.h>
#include <cpl/Mathext.h>
#include "SignalizerDesign.h"
#include <cpl/rendering/OpenGLRasterizers.h>
#include <cpl/simd.h>

namespace Signalizer
{
	// swapping the right channel might give an more intuitive view

	//static const char * ChannelDescriptions[] = { "+R", "+L", "-R", "-L" };
	
	static const char * ChannelDescriptions[] = { "+L", "+R", "-L", "-R" };
	static std::vector<std::string> OperationalModeNames = {"Lissajous", "Polar"};
	
	enum class OperationalModes
	{
		Lissajous,
		Polar
		
	};
	
	enum Textures
	{
		LPlus,
		RPlus,
		LMinus,
		RMinus
	};


	template<typename T, std::size_t size>
	class LookupTable
	{
	public:
		typedef T Ty;
		static const std::size_t tableSize = size;

		inline Ty linearLookup(Ty dx) const noexcept
		{
			Ty scaled = dx * tableSize;
			std::size_t x1 = std::size_t(scaled);
			std::size_t x2 = x1 + 1;
			Ty fraction = scaled - x1;

			return table[x1] * (Ty(1) - fraction) + table[x2] * fraction;
		}

		Ty table[tableSize + 1];
	};
	template<typename T, std::size_t size>
		class QuarterCircleLut : public LookupTable<T, size>
		{
		public:
			QuarterCircleLut()
			{
				double increase = 1.0 / (size - 1);
				for (int i = 0; i < size; ++i)
				{
					
					// describe frist left upper part of circle
					this->table[i] = (T)std::sin(std::acos(1.0 - increase * i));
				}
				this->table[this->tableSize] = (T)1;
			}


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
				section->addControl(&krotation, 0);
				section->addControl(&kwindow, 1);
				section->addControl(&kgain, 0);
				section->addControl(&kenvelopeFollow, 1);
				section->addControl(&kopMode, 0);
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
				section->addControl(&kdiagnostics, 3);
				page->addSection(section, "Options");
			}
			if (auto section = new Signalizer::CContentPage::MatrixSection())
			{
				section->addControl(&kdrawingColour, 0);
				section->addControl(&kgraphColour, 0);
				section->addControl(&kbackgroundColour, 0);
				section->addControl(&kskeletonColour, 0);
				section->addControl(&kprimitiveSize, 1);
				page->addSection(section, "Look");
			}
		}
		return std::unique_ptr<juce::Component>(content);

	}

	CVectorScope::CVectorScope(cpl::AudioBuffer & data)
	:
		audioStream(data),
		kwindow("Window size", cpl::CKnobSlider::ControlType::ms),
		krotation("Wave Z-rotation"),
		kgain("Input gain"),
		kgraphColour("Graph colour"),
		kbackgroundColour("Background colour"),
		kdrawingColour("Drawing colour"),
		kskeletonColour("Skeleton colour"),
		kprimitiveSize("Primitive size"),
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
		krotation.bAddFormatter(this);
		kprimitiveSize.bAddFormatter(this);
		kenvelopeFollow.bAddPassiveChangeListener(this);
		kopMode.bAddPassiveChangeListener(this);
		// buttons
		kantiAlias.bSetTitle("Antialias");
		kantiAlias.setToggleable(true);
		kfadeOld.bSetTitle("Fade older points");
		kfadeOld.setToggleable(true);
		kdrawLines.bSetTitle("Interconnect samples");
		kdrawLines.setToggleable(true);
		kdiagnostics.bSetTitle("Diagnostics");
		kdiagnostics.setToggleable(true);
		kenvelopeFollow.bSetTitle("Auto gain");
		kenvelopeFollow.setToggleable(true);
		
		// design
		kopMode.setValues(OperationalModeNames);
		kopMode.bSetTitle("Operational mode");
		
		
		// descriptions.
		kwindow.bSetDescription("The size of the displayed time window.");
		kgain.bSetDescription("How much the input (x,y) is scaled (or the input gain)" \
							 " - additional transform that only affects the waveform, and not the graph");
		krotation.bSetDescription("The amount of degrees to rotate the waveform around the origin (z-rotation)"\
			" - additional transform that only affects the waveform, and not the graph.");
		kantiAlias.bSetDescription("Antialiases rendering (if set - see global settings for amount). May slow down rendering.");
		kfadeOld.bSetDescription("If set, gradually older samples will be faded linearly.");
		kdrawLines.bSetDescription("If set, interconnect samples linearly.");
		kdrawingColour.bSetDescription("The main colour to paint with.");
		kgraphColour.bSetDescription("The colour of the graph.");
		kbackgroundColour.bSetDescription("The background colour of the view.");
		kdiagnostics.bSetDescription("Toggle diagnostic information in top-left corner.");
		kskeletonColour.bSetDescription("The colour of the box skeleton indicating the OpenGL camera clip box.");
		kprimitiveSize.bSetDescription("The size of the rendered primitives (eg. lines or points).");
		kenvelopeFollow.bSetDescription("Monitors the audio stream and automatically scales the input gain such that it approaches unity intensity");
		kopMode.bSetDescription("Changes the presentation of the data - Lissajous is the classic XY mode on oscilloscopes, while the polar mode is a wrapped circle of the former.");
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
		archive << kdiagnostics;
		archive << kdrawLines;
		archive << kgraphColour;
		archive << kbackgroundColour;
		archive << kdrawingColour;
		archive << ktransform.getTransform3D();
		archive << kskeletonColour;
		archive << kprimitiveSize;
	}

	void CVectorScope::load(cpl::CSerializer::Builder & builder, long long int version)
	{
		builder >> kwindow;
		builder >> kgain;
		builder >> krotation;
		builder >> kantiAlias;
		builder >> kfadeOld;
		builder >> kdiagnostics;
		builder >> kdrawLines;
		builder >> kgraphColour;
		builder >> kbackgroundColour;
		builder >> kdrawingColour;
		builder >> ktransform.getTransform3D();
		builder >> kskeletonColour;
		builder >> kprimitiveSize;

		ktransform.syncEditor();

	}

	void CVectorScope::paint(juce::Graphics & g)
	{

		
		auto cStart = cpl::Misc::ClockCounter();
		

		
		// do software rendering
		if(!isOpenGL())
		{
			g.fillAll(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f));
			g.setColour(kbackgroundColour.getControlColourAsColour().withAlpha(1.0f).contrasting());
			g.drawText("Enable OpenGL in settings to use the vectorscope", getLocalBounds(), juce::Justification::centred);
			
			// post fps
			auto tickNow = juce::Time::getHighResolutionTicks();
			avgFps.setNext(tickNow - lastFrameTick);
			lastFrameTick = tickNow;
			
		}
		
		if (kdiagnostics.bGetValue() > 0.5)
		{
			auto fps = 1.0 / (avgFps.getAverage() / juce::Time::getHighResolutionTicksPerSecond());

			auto totalCycles = renderCycles + cpl::Misc::ClockCounter() - cStart;
			double cpuTime = (double(totalCycles) / (processorSpeed * 1000 * 1000) * 100) * fps;
			g.setColour(juce::Colours::blue);
			sprintf(textbuf.get(), "%dx%d: %.1f fps - %.1f%% cpu", 
				getWidth(), getHeight(), fps, cpuTime);
			g.drawSingleLineText(textbuf.get(), 10, 20);
			
		}
	}

	void CVectorScope::initOpenGL()
	{
		const int imageSize = 128;
		const float fontToPixelScale = 90 / 64.0f;
		
		textures.clear();

		for (int i = 0; i < ArraySize(ChannelDescriptions); ++i)
		{

			juce::Image letter(juce::Image::ARGB, imageSize, imageSize, true);
			{
				juce::Graphics g(letter);
				g.fillAll(juce::Colours::transparentBlack);
				g.setColour(Colours::white);
				g.setFont(letter.getHeight() * fontToPixelScale * 0.5);
				g.drawText(ChannelDescriptions[i], letter.getBounds().toFloat(), juce::Justification::centred, false);
			}
			textures.push_back(std::unique_ptr<juce::OpenGLTexture>((new juce::OpenGLTexture())));
			textures[i]->loadImage(letter);
		}
	}
	void CVectorScope::closeOpenGL()
	{
		textures.clear();
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
		const float psize = kprimitiveSize.bGetValue() * 10;
		const float gain = getGain();
		const bool isPolar = kopMode.getZeroBasedSelIndex() == (int)OperationalModes::Polar;

		float sleft(0.0f), sright(0.0f);

		cpl::AudioBuffer * buffer;

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
		const std::size_t numSamples = (*buffer)[0].size - 1;

		cpl::OpenGLEngine::COpenGLStack openGLStack;
		// set up openGL
		openGLStack.setBlender(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
		openGLStack.loadIdentityMatrix();


		openGLStack.applyTransform3D(ktransform.getTransform3D());
		antiAlias ? openGLStack.enable(GL_MULTISAMPLE) : openGLStack.disable(GL_MULTISAMPLE);

		cpl::CAudioBuffer::CChannelBuffer * leftChannel, * rightChannel;
		leftChannel = &((*buffer)[0]);
		rightChannel = &((*buffer)[1]);

		openGLStack.setLineSize(psize);
		openGLStack.setPointSize(psize);

		// add a stack scope for transformations.
		if (isPolar)
		{
			using namespace cpl::simd;
			cpl::OpenGLEngine::MatrixModification matrixMod;
			// fixate the rotation to fit.
			//matrixMod.rotate(135, 0, 0, 1);
			// and apply the gain:
			matrixMod.scale(gain, gain, 1);

			cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);

			drawer.addColour(kdrawingColour.getControlColourAsColour());

			static const v4sf vPiHalf = set1<v4sf>(M_PI / 2);
			static const v4sf vSignBit = set1<v4sf>(-0.0f);
			static const v4sf vHalf = set1<v4sf>(0.5f);
			static const v4sf vOne = set1<v4sf>(1.0f);
			static const v4sf vZero = set1<v4sf>(0.0f);
			static const v4sf vCosine = set1<v4sf>((float)std::cos(M_PI * 135.0 / 180));
			static const v4sf vSine = set1<v4sf>((float)std::sin(M_PI * 135.0 / 180));

			float sampleFade = 1.0 / numSamples;

			suitable_container<v4sf> outX, outY;

			// iterate front of buffer, then back.
			// this feels awkward, but that's kinda how circular 
			// buffers work.
			std::size_t depthCounter = 0;
			for (int section = 0; section < 2; ++section)
			{
				auto const lit = leftChannel->getIterator<4>();
				auto const rit = rightChannel->getIterator<4>();

				auto const sectionSamples = lit.sizes[section];


				for (std::size_t i = 0; i < sectionSamples; i += elements_of<v4sf>::value)
				{
					const v4sf vLeft = loadu<v4sf>(lit.getIndex(section) + i);
					const v4sf vRight = loadu<v4sf>(rit.getIndex(section) + i);




					// rotate our view manually (it needs to be fixed on y axis for the next math to work)
					const v4sf vY = vLeft * vCosine - vRight * vSine;
					const v4sf vX = vLeft * vSine + vRight * vCosine;

					
					
					// the polar display does a max on the y-coordinate (so everything stays over the line)
					// and swaps the sign of the x-coordinate if y crosses zero.
					auto const vSignMask = vand((vY < vZero), vSignBit);
					outX = vxor(vX, vSignMask); // x = y < 0 ? -x : x
					outY = abs(vY); // y = |y|

					// draw vertices
					drawer.addVertex(outX[0], outY[0], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[1], outY[1], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[2], outY[2], depthCounter++ * sampleFade - 1);
					drawer.addVertex(outX[3], outY[3], depthCounter++ * sampleFade - 1);
				}
				
				/*for (std::size_t i = 0; i < sectionSamples; i += elements_of<v4sf>::value)
				 {
				 const v4sf vLeft = loadu<v4sf>(lit.getIndex(section) + i);
				 const v4sf vRight = loadu<v4sf>(rit.getIndex(section) + i);
				 
				 
				 
				 
				 // rotate our view manually (it needs to be fixed on y axis for the next math to work)
				 const v4sf vY = vLeft * vCosine - vRight * vSine;
				 const v4sf vX = vLeft * vSine + vRight * vCosine;
				 
				 //cart to polar:
				// abs(x,y) * e^i*atan2f(x, y)
				
					 
					 
				 
				 // the polar display does a max on the y-coordinate (so everything stays over the line)
				 // and swaps the sign of the x-coordinate if y crosses zero.
				 auto const vSignMask = vand((vY < vZero), vSignBit);
				 outX = cpl::sse::sin_ps(vX); // x = y < 0 ? -x : x
				outY = cpl::sse::sin_ps(vY); // y = |y|
				 
				 // draw vertices
				 drawer.addVertex(outX[0], outY[0], depthCounter++ * sampleFade - 1);
				 drawer.addVertex(outX[1], outY[1], depthCounter++ * sampleFade - 1);
				 drawer.addVertex(outX[2], outY[2], depthCounter++ * sampleFade - 1);
				 drawer.addVertex(outX[3], outY[3], depthCounter++ * sampleFade - 1);
				 }*/
				
			}


		}
		else
		{
			cpl::OpenGLEngine::MatrixModification matrixMod;
			// apply the custom rotation to the waveform
			matrixMod.rotate(krotation.bGetValue() * 360, 0, 0, 1);
			// and apply the gain:
			matrixMod.scale(gain, gain, 1);

			float sampleFade = 1.0 / numSamples;

			if (!fadeHistory)
			{
				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);

				drawer.addColour(kdrawingColour.getControlColourAsColour());
				for (std::size_t i = 0; i < numSamples; ++i)
				{
					// use glDrawArrays instead here
					sleft = leftChannel->singleCheckAccess(i);
					sright = rightChannel->singleCheckAccess(i);

					drawer.addVertex(sright, sleft, i * sampleFade - 1);
				}
			}
			else
			{
				cpl::OpenGLEngine::PrimitiveDrawer<1024> drawer(openGLStack, fillPath ? GL_LINE_STRIP : GL_POINTS);
				
				auto jcolour = kdrawingColour.getControlColourAsColour();
				float red = jcolour.getFloatRed(), blue = jcolour.getFloatBlue(),
					green = jcolour.getFloatGreen(), alpha = jcolour.getFloatGreen();

				float fade = 0;

				for (std::size_t i = 0; i < numSamples; ++i)
				{
					sleft = leftChannel->singleCheckAccess(i);
					sright = rightChannel->singleCheckAccess(i);

					fade = i * sampleFade;

					drawer.addColour(fade * red, fade * green, fade * blue, alpha);
					drawer.addVertex(sright, sleft, i * sampleFade - 1);
				
				}
			}
		}

		openGLStack.setLineSize(2.0f);

		// draw skeleton graph
		if (!isPolar)
		{
			cpl::OpenGLEngine::PrimitiveDrawer<128> drawer(openGLStack, GL_LINES);
			drawer.addColour(kskeletonColour.getControlColourAsColour());
			int nlines = 14;
			auto rel = 1.0f / nlines;

			// front vertival
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(i * rel * 2 - 1, -1.f, 0.0f);
				drawer.addVertex(i * rel * 2 - 1, 1.f, 0.0f);
			}
			// front horizontal
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(-1.0f, i * rel * 2 - 1, 0.0f);
				drawer.addVertex(1.0f, i * rel * 2 - 1, 0.0f);
			}
			// back vertical
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(i * rel * 2 - 1, -1.f, -1.0f);
				drawer.addVertex(i * rel * 2 - 1, 1.f, -1.0f);
			}
			// back horizontal
			for (int i = 0; i <= nlines; ++i)
			{
				drawer.addVertex(-1.0f, i * rel * 2 - 1, -1.0f);
				drawer.addVertex(1.0f, i * rel * 2 - 1, -1.0f);
			}
		}
		else
		{
			// draw two half circles
			auto lut = circleData.get();
			
			int numInt = lut->tableSize;
			float advance = 1.0f / (numInt - 1);
			{
				cpl::OpenGLEngine::PrimitiveDrawer<512> drawer(openGLStack, GL_LINES);
				drawer.addColour(kskeletonColour.getControlColourAsColour());

				float oldY = 0.0f;
				for (int i = 1; i < numInt; ++i)
				{
					auto fraction = advance * i;
					auto yCoordinate = lut->linearLookup(fraction);
					auto leftX = -1.0f + fraction;
					auto rightX = 1.0f - fraction;
					// left part
					drawer.addVertex(leftX - advance, oldY, 0);
					drawer.addVertex(leftX, yCoordinate, 0);
					drawer.addVertex(leftX - advance, oldY, -1);
					drawer.addVertex(leftX, yCoordinate, -1);
					
					// right part
					drawer.addVertex(rightX + advance, oldY, 0);
					drawer.addVertex(rightX, yCoordinate, 0);
					drawer.addVertex(rightX + advance, oldY, -1);
					drawer.addVertex(rightX, yCoordinate, -1);
					//drawer.addVertex(1 - fraction, yCoordinate, 0);
					
					oldY = yCoordinate;
				}
				
				// add front and back horizontal lines.
				drawer.addVertex(-1.0f, 0.0f, 0.0f);
				drawer.addVertex(1.0f, 0.0f, 0.0f);
				drawer.addVertex(-1.0f, 0.0f, -1.0f);
				drawer.addVertex(1.0f, 0.0f, -1.0f);
				
				// add critical diagonal phase lines.
				const float quarterPISinCos = 0.707106781186547f;
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(quarterPISinCos, quarterPISinCos, 0.0f);
				drawer.addVertex(0.0f, 0.0f, 0.0f);
				drawer.addVertex(-quarterPISinCos, quarterPISinCos, 0.0f);
				
				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(quarterPISinCos, quarterPISinCos, -1.0f);
				drawer.addVertex(0.0f, 0.0f, -1.0f);
				drawer.addVertex(-quarterPISinCos, quarterPISinCos, -1.0f);;
			}
		}

		// Draw basic graph

		{
			cpl::OpenGLEngine::PrimitiveDrawer<12> drawer(openGLStack, GL_LINES);
			drawer.addColour(kgraphColour.getControlColourAsColour());
			// front x, y axii
			drawer.addVertex(-1.0f, 0.0f, 0.0f);
			drawer.addVertex(1.0f, 0.0f, 0.0f);
			drawer.addVertex(0.0f, 1.0f, 0.0f);
			drawer.addVertex(0.0f, isPolar ? 0.0f : -1.0f, 0.0f);
			

		}


		// draw channel rotations letters.
		{
			auto rotation = -krotation.bGetValue() * 2 * M_PI;
			const float heightToWidthFactor = float(getHeight()) / getWidth();
			using namespace cpl::simd;

			cpl::OpenGLEngine::MatrixModification m;
			// this undoes the text squashing due to variable aspect ratios.
			m.scale(heightToWidthFactor, 1.0f, 1.0f);

			// calculate coordinates using sin/cos simd pairs.
			suitable_container<v4sf> phases, xcoords, ycoords;

			// set phases (we rotate L/R etc. text around in a circle)
			phases[0] = rotation;
			phases[1] = M_PI * 0.5f + rotation;
			phases[2] = M_PI + rotation;
			phases[3] = M_PI * 1.5f + rotation;

			const float circleScaleFactor = 1.1f;

			// some registers
			v4sf 
				vsines, 
				vcosines, 
				vscale = set1<v4sf>(circleScaleFactor), 
				vadd = set1<v4sf>(1.0f - circleScaleFactor),
				vheightToWidthFactor = set1<v4sf>(1.0f / heightToWidthFactor);

			// do 8 trig functions in one go!
			cpl::sse::sincos_ps(phases, &vsines, &vcosines);

			// place the circle just outside the graph, and offset it.
			xcoords = vsines * vscale * vheightToWidthFactor + vadd;
			ycoords = vcosines * vscale + vadd;

			auto jcolour = kgraphColour.getControlColourAsColour();
			// render texture text at coordinates.
			for (int i = 0; i < 4; ++i)
			{
				cpl::OpenGLEngine::ImageDrawer text(openGLStack, *textures[i]);
				text.setColour(jcolour);
				text.drawAt({ xcoords[i], ycoords[i], 0.1f, 0.1f });
			}
		}

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
		}
		else
		{
			matrix.position.x += deltaDifference.x / 500.f;
			matrix.position.y += factor * -deltaDifference.y / 500.f;
		}
		ktransform.syncEditor();
		
		
		lastMousePos = event.position;
	}
	void CVectorScope::mouseUp(const MouseEvent& event)
	{

	}
	void CVectorScope::mouseDown(const MouseEvent& event)
	{
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
		else if (ctrl == &krotation)
		{
			sprintf(buf, "%.2f degs", value * 360);
			buffer = buf;
			return true;
		}
		else if (ctrl == &kprimitiveSize)
		{
			sprintf(buf, "%.2f pts", value * 10);
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
		else if (ctrl == &krotation)
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
		else if (ctrl == &kprimitiveSize)
		{
			double newValue;
			char * endPtr = nullptr;
			newValue = strtod(buffer.c_str(), &endPtr);
			if (endPtr > 0)
			{
				value = cpl::Math::confineTo(newValue / 10.0, 0.0, 1.0);
				return true;
			}
		}
		return false;
	}
	void CVectorScope::onObjectDestruction(const cpl::CBaseControl::ObjectProxy & destroyedObject)
	{
		// hmmm.....
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
	}


};
