#include <cmath>
#include <fstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <ppl.h>
#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"
#include "timer.hpp"

class FlockEvoApp {

public:
  FlockEvoApp();
  void run();

private:
	void simulateGeneration(Predator &predator, std::vector<Prey> tempPrey);
	void createNewGeneration();

	std::vector<Predator> predatorAnimats;
	std::vector<std::vector<Prey>> preyAnimats;
	std::vector<std::vector<Prey>> preyAnimats2;
	
	int generationCount;
	int stepCount;
	int currentPredator;
	int totalHuntCount;

	std::ofstream CSVFile, CSVFilePrey;
};

FlockEvoApp::FlockEvoApp() {
	rnd_seed((unsigned)time(NULL));
	srand(time(NULL));
}

void FlockEvoApp::run() {

  // for all iterations
  for (int iterationCount = 0; iterationCount < AppSettings::noOfIterations; iterationCount++) {
    // init
    generationCount = 0;
    stepCount = 0;
    currentPredator = 0;
	totalHuntCount = 0;

    // log file
    char fileName[64];
    sprintf_s(fileName, 64, "log%d", iterationCount + 1);
    char extension[8] = ".csv";
    strcat_s(fileName, 64, extension);
    CSVFile = std::ofstream(fileName);

	// log file prey
	sprintf_s(fileName, 64, "logPrey%d", iterationCount + 1);
	strcat_s(fileName, 64, extension);
	CSVFilePrey = std::ofstream(fileName);

    // create predators
    predatorAnimats.clear();
    for(int i = 0; i < AppSettings::noOfPredatorAnimats; i++) {
      predatorAnimats.push_back(Predator(i + 1));
    }

	// create prey flocks (save them to preyAnimats)
	preyAnimats.clear();
	preyAnimats2.clear();
	for (int j = 0; j < 500; j++) {
		std::vector<Prey> tempPrey;
		for (int i = 0; i < AppSettings::noOfPreyAnimats; i++) {
			tempPrey.push_back(Prey(i));
		}
		preyAnimats.push_back(tempPrey);
	}

    for (generationCount = 0; generationCount < AppSettings::noOfGenerations; generationCount++) {
      
		// v prey in predAnimats bodo po tem klicu GA rezultati
		if (generationCount != 0) createNewGeneration();
		timer tg;

			// repeat noOfFlocks times
			for (int flockCount = 0; flockCount < AppSettings::noOfFlocks; flockCount++) {

				// premade flocks based on evolution results
				std::vector<Prey> tempPrey = preyAnimats[generationCount * 5 + flockCount];

				std::for_each(predatorAnimats.begin(), predatorAnimats.end(), [this, &tempPrey](Predator& pred){
					simulateGeneration(pred, tempPrey);
				});
			}

			char buffer[128];
			int deadCount = 0;
			#if (PREY_EVOLUTION == 1)
			int count = 0;
				for (std::vector<Prey> prey : preyAnimats2)
				{
					count++;
					for (int k = 0; k < 100; k++) {
						if (prey[k].isDead) deadCount++;
					}

					CSVFilePrey << iterationCount + 1 << ';';
					CSVFilePrey << generationCount + 1 << ';';
					CSVFilePrey << count << ';';
					CSVFilePrey << prey[0].wb << ';';
					CSVFilePrey << prey[0].wr << ';';
					CSVFilePrey << prey[0].we << ';';
					CSVFilePrey << prey[0].cn0 << ';';
					CSVFilePrey << prey[0].cn1 << ';';
					CSVFilePrey << prey[0].cn2 << ';';
					CSVFilePrey << prey[0].cn3 << ';';
					CSVFilePrey << prey[0].cr0 << ';';
					CSVFilePrey << prey[0].cr1 << ';';
					CSVFilePrey << prey[0].selfishTime << ';';
					CSVFilePrey << prey[0].selfishTime2 << ';';
					CSVFilePrey << prey[0].selfishEscapeDistance << ';';
					CSVFilePrey << prey[0].selfishProbability << ';';
					CSVFilePrey << prey[0].dot_pow << ';';
					CSVFilePrey << prey[0].dotTreshold << ';';
					CSVFilePrey << prey[0].v_b << ';';
					CSVFilePrey << prey[0].v_r << ';';
					CSVFilePrey << std::endl;
				}
				
				//sprintf_s(buffer, "%d %d %d %g s\n", iterationCount + 1, generationCount + 1, tg.elapsed());
				//std::cout << buffer;
			#endif

				int reportHuntCount = 0;
				for (Predator const& predator : predatorAnimats)
				{
					totalHuntCount += predator.huntCount;

					//int reportHuntCount = (predator.huntCount - AppSettings::initialHuntCount) / AppSettings::huntFactor;
					reportHuntCount += predator.huntCount;
					CSVFile << iterationCount + 1 << ';';
					CSVFile << generationCount + 1 << ';';
					CSVFile << predator.id << ';';
					CSVFile << predator.nearestWeight << ';';
					CSVFile << predator.peripheralWeight << ';';
					CSVFile << predator.centreWeight << ';';
					CSVFile << predator.distanceForAcceleration << ';';
					CSVFile << predator.velocityMultiplier << ';';
					CSVFile << predator.attackPeriod << ';';
					CSVFile << predator.huntCount;
					CSVFile << std::endl;
				
					//sprintf_s(buffer, "%d;%d;%d;%.2f;%.2f;%.2f;%d\n", iterationCount + 1, generationCount + 1, predator.id, predator.nearestWeight, predator.peripheralWeight, predator.centreWeight, reportHuntCount);
					//CSVFile << buffer;
				}
  
				sprintf_s(buffer, "%d %d %d %d %g s\n", iterationCount + 1, generationCount + 1, reportHuntCount, deadCount, tg.elapsed());
				std::cout << buffer;

				CSVFile.close();
				CSVFilePrey.close();
			
    }
  }
}

// search pairwise neighbours
void searchNeighbours(std::vector<Prey>& prey)
{
  const float maxSize = std::max(AppSettings::cohesionSize, AppSettings::peripheralitySize);
  const float maxSize2 = maxSize * maxSize;
  auto end = prey.end();
  for (auto a = prey.begin(); a != end; ++a)
  {
    if (!a->isDead)
    {
      for (auto b = a + 1; b != end; ++b)
      {
        if (!b->isDead)
        {
          const glm::vec2 ofs = b->position - a->position;
          const float dist2 = glm::length2(ofs);
          if ((dist2 > 0.0f) && (dist2 < maxSize2))
          {
						float dist = std::sqrt(dist2) - (2 * AppSettings::preySize);
						if (dist < .0f)
							dist = .0f;
						if (dist < maxSize)
						{
							a->neighbours.emplace_back(ofs, b->getVelocity(), dist);
							b->neighbours.emplace_back(-ofs, a->getVelocity(), dist);
						}
          }
        }
      }
    }
  }
}

void FlockEvoApp::simulateGeneration(Predator &predator, std::vector<Prey> tempPrey)
{
  for (int i = 0; i < AppSettings::noOfSteps; i++)
  {
    // calculate
    // prey
    searchNeighbours(tempPrey);
    for(std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
    {
      if (!p->isDead) p->calculate(predator);
    }
    // predator
    predator.calculate(tempPrey);

    // update
    // prey
    for(std::vector<Prey>::iterator p = tempPrey.begin(); p != tempPrey.end(); ++p)
    {
      if (!p->isDead)
        p->update(predator, tempPrey);
    }
    // predator
    predator.update(tempPrey);
  }

	preyAnimats2.push_back(tempPrey);
	predator.reset();
}

void FlockEvoApp::createNewGeneration() {

  // predator
  float x = AppSettings::worldCentre.x;
  float y = AppSettings::worldCentre.y + (AppSettings::preySize * 200.0f);

  float vX = .0f;
  float vY = -1.0f;

  //genetic stuff
  std::vector<Predator> oldPredators = predatorAnimats;
  predatorAnimats.clear();
  for(int i = 0; i < AppSettings::noOfPredatorAnimats; i++)
  {
    // first select 2 random parents
    int randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
    Predator firstParent = oldPredators.at(randomIndex);
    randomIndex = randomInt(0, AppSettings::noOfPredatorAnimats - 1);
    Predator secondParent = oldPredators.at(randomIndex);

    // get parents
    if (totalHuntCount > 1)
    {
      int currentHuntCount = 0;
      int randomHuntCount = randomInt(1, totalHuntCount);
      // first parent
      for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++)
      {
        firstParent = oldPredators.at(j);
        currentHuntCount += firstParent.huntCount;
        if (currentHuntCount >= randomHuntCount)
          break;
      }

      currentHuntCount = 0;
      randomHuntCount = randomInt(1, totalHuntCount);
      // second parent
      for (int j = 0; j < AppSettings::noOfPredatorAnimats; j++)
      {
        secondParent = oldPredators.at(j);
        currentHuntCount += secondParent.huntCount;
        if (currentHuntCount >= randomHuntCount)
          break;
      }
    }

	// coinflip crossover
	float aP; 
	float vM;
	float dA;

	float random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) aP = firstParent.attackPeriod;
	else aP = secondParent.attackPeriod;

	random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) vM = firstParent.velocityMultiplier;
	else vM = secondParent.velocityMultiplier;

	random = randomFloat(.0f, 1.0f);
	if (random < 0.5f) dA = firstParent.distanceForAcceleration;
	else dA = secondParent.distanceForAcceleration;

	// mutation

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) aP = aP + (aP * AppSettings::mutationFactor);
		else aP = aP - (aP * AppSettings::mutationFactor);
	}
	if (aP < 0) aP = 0;
	if (aP > AppSettings::noOfSteps) aP = AppSettings::noOfSteps;

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) vM += AppSettings::mutationFactor;
		else vM -= AppSettings::mutationFactor;
	}
	if (vM < 1.0f) vM = 1.0f;
	if (vM > 3.0f) vM = 3.0f;

	random = randomFloat(.0f, 1.0f);
	if (random < AppSettings::mutationRate){
		float flip = randomFloat(.0f, 1.0f);
		if (flip > 0.5f) dA = dA + (dA * AppSettings::mutationFactor);
		else dA = dA - (dA * AppSettings::mutationFactor);;
	}
	if (dA < 0) dA = 0;
	if (dA > AppSettings::huntSize) dA = AppSettings::huntSize;

    // coinflip crossover
    float nW;
    float pW;
    float cW;
    random = randomFloat(.0f, 1.0f);
    if (random < 0.5f)
      nW = firstParent.nearestWeight;
    else
      nW = secondParent.nearestWeight;
    random = randomFloat(.0f, 1.0f);
    if (random < 0.5f)
      pW = firstParent.peripheralWeight;
    else
      pW = secondParent.peripheralWeight;
    random = randomFloat(.0f, 1.0f);
    if (random < 0.5f)
      cW = firstParent.centreWeight;
    else
      cW = secondParent.centreWeight;

    // mutation
    random = randomFloat(.0f, 1.0f);
    if (random < AppSettings::mutationRate)
    {
      float flip = randomFloat(.0f, 1.0f);
      if (flip > 0.5f)
        nW += AppSettings::mutationFactor;
      else
        nW -= AppSettings::mutationFactor;
    }
    random = randomFloat(.0f, 1.0f);
    if (random < AppSettings::mutationRate)
    {
      float flip = randomFloat(.0f, 1.0f);
      if (flip > 0.50f)
        pW += AppSettings::mutationFactor;
      else
        pW -= AppSettings::mutationFactor;
    }
    random = randomFloat(.0f, 1.0f);
    if (random < AppSettings::mutationRate)
    {
      float flip = randomFloat(.0f, 1.0f);
      if (flip > 0.5f)
        cW += AppSettings::mutationFactor;
      else
        cW -= AppSettings::mutationFactor;
    }

    if (nW < 0)
      nW = 0;
    if (pW < 0)
      pW = 0;
    if (cW < 0)
      cW = 0;

    // weigh
    float sum = .0f;
    sum += nW;
    sum += pW;
    sum += cW;
    nW /= sum;
    pW /= sum;
    cW = 1 - (nW + pW);

    predatorAnimats.push_back(Predator(i + 1, nW, pW, cW, dA, vM, aP));
  }

  currentPredator = 0;
  totalHuntCount = 0;

#if (PREY_EVOLUTION == 1)

  std::vector<std::vector<Prey>> oldPrey = preyAnimats2;
  preyAnimats2.clear();

  float totalPreyFitness = 0.0f;
  std::vector<int> noAlive;
  std::vector<float> avgEnergy;
  float maxEnergy = 0;
  std::vector<float> preyFitness;

  for (int i = 0; i < oldPrey.size(); i++) {
	  float aliveCount = 0;
	  float totalEnergy = 0;
	  for (int j = 0; j < AppSettings::noOfPreyAnimats; j++) {
		  if (!oldPrey[i][j].isDead) {
			  aliveCount++;
			  totalEnergy += oldPrey[i][j].energy;
		  }
	  }
	  noAlive.push_back(aliveCount);
	  avgEnergy.push_back(totalEnergy / aliveCount);
	  if (totalEnergy / aliveCount > maxEnergy) maxEnergy = totalEnergy / aliveCount;
  }

  // get fitness for each flock
  for (int i = 0; i < oldPrey.size(); i++) {
	  preyFitness.push_back(noAlive[i] / AppSettings::noOfPreyAnimats + avgEnergy[i] / maxEnergy);
	  totalPreyFitness += preyFitness[i];
  }

  for (int i = 0; i < 500; i++) {
	  
	  int randomIndex = randomInt(0, oldPrey.size() - 1);
	  std::vector<Prey> firstParent = oldPrey.at(randomIndex);

	  randomIndex = randomInt(0, oldPrey.size() - 1);
	  std::vector<Prey> secondParent = oldPrey.at(randomIndex);

	  if (totalPreyFitness > 0) {

		  float currentFitness = 0;
		  float randomFitness = randomFloat(0.f, totalPreyFitness);

		  for (int j = 0; j < oldPrey.size(); j++) {
			  firstParent = oldPrey.at(j);
			  currentFitness += preyFitness[j];
			  if (currentFitness > randomFitness) break;
		  }

		  currentFitness = 0;
		  randomFitness = randomFloat(0.f, totalPreyFitness);

		  for (int j = 0; j < oldPrey.size(); j++) {
			  secondParent = oldPrey.at(j);
			  currentFitness += preyFitness[j];
			  if (currentFitness >= randomFitness) break;
		  }
	  }

	  float wb1, wr1, we1;
	  float cn01, cn11, cn21, cn31;
	  float cr01, cr11;
	  int selfishTime1, selfishTime21;
	  float selfishEscapeDistance1, selfishProbability1;
	  float dot_pow1, dotTreshold1;
	  float v_b1, v_r1;

	  // coinflip crossover

	  float random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) wb1 = firstParent[0].wb;
	  else wb1 = secondParent[0].wb;
	  
	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) wr1 = firstParent[0].wr;
	  else wr1 = secondParent[0].wr;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) we1 = firstParent[0].we;
	  else we1 = secondParent[0].we;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cn01 = firstParent[0].cn0;
	  else cn01 = secondParent[0].cn0;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cn11 = firstParent[0].cn1;
	  else cn11 = secondParent[0].cn1;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cn21 = firstParent[0].cn2;
	  else cn21 = secondParent[0].cn2;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cn31 = firstParent[0].cn3;
	  else cn31 = secondParent[0].cn3;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cr01 = firstParent[0].cr0;
	  else cr01 = secondParent[0].cr0;
	  
	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) cr11 = firstParent[0].cr1;
	  else cr11 = secondParent[0].cr1;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) selfishTime1 = firstParent[0].selfishTime;
	  else selfishTime1 = secondParent[0].selfishTime;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) selfishTime21 = firstParent[0].selfishTime2;
	  else selfishTime21 = secondParent[0].selfishTime2;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) selfishEscapeDistance1 = firstParent[0].selfishEscapeDistance;
	  else selfishEscapeDistance1 = secondParent[0].selfishEscapeDistance;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) selfishProbability1 = firstParent[0].selfishProbability;
	  else selfishProbability1 = secondParent[0].selfishProbability;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) dot_pow1 = firstParent[0].dot_pow;
	  else dot_pow1 = secondParent[0].dot_pow;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) dotTreshold1 = firstParent[0].dotTreshold;
	  else dotTreshold1 = secondParent[0].dotTreshold;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) v_r1 = firstParent[0].v_r;
	  else v_r1 = secondParent[0].v_r;

	  random = randomFloat(0.0f, 1.0f);
	  if (random < 0.5f) v_b1 = firstParent[0].v_b;
	  else v_b1 = secondParent[0].v_b;
	  cr11 = glm::clamp(cr11, 0.0f, 50.0f);

	  // mutation
	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) wb1 += AppSettings::mutationFactor;
		  else wb1 -= AppSettings::mutationFactor;
	  }
	  wb1 = glm::clamp(wb1, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) wr1 += AppSettings::mutationFactor;
		  else wr1 -= AppSettings::mutationFactor;
	  }
	  wr1 = glm::clamp(wr1, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) we1 += AppSettings::mutationFactor;
		  else we1 -= AppSettings::mutationFactor;
	  }
	  we1 = glm::clamp(we1, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cn01 += AppSettings::mutationFactor;
		  else cn01 -= AppSettings::mutationFactor;
	  }
	  cn01 = glm::clamp(cn01, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cn11 += AppSettings::mutationFactor;
		  else cn11 -= AppSettings::mutationFactor;
	  }
	  cn11 = glm::clamp(cn11, 0.0f, 50.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cn21 += AppSettings::mutationFactor;
		  else cn21 -= AppSettings::mutationFactor;
	  }
	  cn21 = glm::clamp(cn21, AppSettings::minPreyVelocity, AppSettings::maxPreyVelocity);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cn31 += AppSettings::mutationFactor;
		  else cn31 -= AppSettings::mutationFactor;
	  }
	  cn31 = glm::clamp(cn31, 0.0f, 50.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cr01 += AppSettings::mutationFactor;
		  else cr01 -= AppSettings::mutationFactor;
	  }
	  cr01 = glm::clamp(cr01, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) cr11 += AppSettings::mutationFactor;
		  else cr11 -= AppSettings::mutationFactor;
	  }
	  cr11 = glm::clamp(cr11, 0.0f, 50.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) selfishTime1 += AppSettings::mutationFactor;
		  else selfishTime1 -= AppSettings::mutationFactor;
	  }
	  selfishTime1 = glm::clamp(selfishTime1, 1, 600);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) selfishTime21 += AppSettings::mutationFactor;
		  else selfishTime21 -= AppSettings::mutationFactor;
	  }
	  selfishTime21 = glm::clamp(selfishTime21, 0, 600);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) selfishEscapeDistance1 += AppSettings::mutationFactor;
		  else selfishEscapeDistance1 -= AppSettings::mutationFactor;
	  }
	  selfishEscapeDistance1 = glm::clamp(selfishEscapeDistance1, 0.0f, AppSettings::huntSize);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) selfishProbability1 += AppSettings::mutationFactor;
		  else selfishProbability1 -= AppSettings::mutationFactor;
	  }
	  selfishProbability1 = glm::clamp(selfishProbability1, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) dot_pow1 += AppSettings::mutationFactor;
		  else dot_pow1 -= AppSettings::mutationFactor;
	  }
	  dot_pow1 = glm::clamp(dot_pow1, 0.0f, 50.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) dotTreshold1 += AppSettings::mutationFactor;
		  else dotTreshold1 -= AppSettings::mutationFactor;
	  }
	  dotTreshold1 = glm::clamp(dotTreshold1, 0.0f, 1.0f);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) v_r1 += AppSettings::mutationFactor;
		  else v_r1 -= AppSettings::mutationFactor;
	  }
	  v_r1 = glm::clamp(v_r1, AppSettings::minPreyVelocity, AppSettings::cruisingPrey);

	  random = randomFloat(0.0f, 1.0f);
	  if (random < AppSettings::mutationRate){
		  float flip = randomFloat(.0f, 1.0f);
		  if (flip > 0.5f) v_b1 += AppSettings::mutationFactor;
		  else v_b1 -= AppSettings::mutationFactor;
	  }
	  v_b1 = glm::clamp(v_b1, AppSettings::cruisingPrey, AppSettings::maxPreyVelocity);

	  std::vector<Prey> tempPrey;
	  for (int j = 0; j < AppSettings::noOfPreyAnimats; j++){
		  tempPrey.push_back(Prey(j, wb1, wr1, we1, cn01, cn11, cn21, cn31, cr01, cr11, selfishTime1, selfishTime21, selfishEscapeDistance1, selfishProbability1, dot_pow1, dotTreshold1, v_b1, v_r1));
	  }
	  preyAnimats.push_back(tempPrey);
	  
  }

#endif
}

int main()
{
  FlockEvoApp app;
  app.run();
}