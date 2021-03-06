#pragma once
#include "glm.hpp"
#include "random.hpp"

namespace AppSettings {

	// best conf 0 = BP
	// best conf 1 = DP
	// best conf 2 = BP

	// 0 - no energy, 1 - oxygen consumption, 2 - zheng
	#define ENERGY 2
	#define PREY_ENERGY 2
	// 0 - predator non confusable, 1 - predator confusable, 2 - predator zheng
	#define CONFUSABILITY 2
	// 0 - no selfish prey escape, 1 - selfish prey escape
	#define SELFISH_ESCAPE 1
	// 0 - no hydrodinamics, 1 - hydrodinamics
	#define HYDRO 0
	// 0 - no evolution parameters in account, 1 - evolution parameters present
	#define PREY_EVOLUTION 0 // 2

	#define PREY_ENERGY_PARAMS 1 // 2 in 3
	#define PRED_ENERGY_PARAMS 1 // 3

	//model
	const int screenWidth = 1920;
	const int screenHeight = 1080;

	const int noOfIterations = 20; // 20
	const int noOfSteps = 600; // 600
	const int noOfGenerations = 500; // 500
	const int noOfFlocks = 5; // 5
	const float mutationRate = .02f;
	const float mutationFactor = .2f;
	const float mutationFactorInt = 1;

	const float worldSize = 100;
	const glm::vec2 worldCentre = glm::vec2((float)screenWidth/2, (float)screenHeight/2);

	// prey settings
	const int noOfPreyAnimats = 100;
	const float preySize = 1.0f;

	const float maxPreyForce = 2.0f;
	const float maxPreyVelocity = preySize * 4.0f;
	const float cruisingPrey = preySize * 1.4f;
	const float minPreyVelocity = preySize * 1.0f;

	const float escapeSize = 200.0f;
	const float separationSize = 10.0f;
	const float alignmentSize = 50.0f;
	const float cohesionSize = 200.0f;
	const float peripheralitySize = 200.0f;

	const float escapeWeight = 5.0f;
	const float separationWeight = 5.0f;
	const float alignmentWeight = .3f;
	const float cohesionWeight = .01f;

	// full blind angle (deg)
	const float preyBlindAngle = 60.0f;
	const float cosPreyBlindThreshold = glm::cos(glm::radians(180.0f - preyBlindAngle / 2.0f));

	// predator settings
	const int noOfPredatorAnimats = 100;
	const float predatorSize = 3.0f;
	const int handlingTime = 30;
	const int initialHuntCount = 1;
	const int huntFactor = 1;

	const float maxPredatorForce = 2.5f;
	const float maxPredatorVelocity = preySize * 5.6f;
	const float cruisingSpeed = predatorSize * 1.4f; // brainbridge
	const float minPredatorVelocity = predatorSize * 1.0f;

	const float huntSize = 800.0f;
	const float confusabilitySize = 50.0f;

	const float catchDistance = 12.0f;
	const float oppositeThreshold = glm::cos(glm::half_pi<float>());
}