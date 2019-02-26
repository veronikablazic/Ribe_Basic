#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"
#include <algorithm>

// tactics
#define NEAREST 0
#define PERIPHERAL 1
#define CENTRE 2

Predator::Predator() {}

Predator::Predator(int animatID) {
  // acceleration
  acceleration = glm::vec2(.0f, .0f);

  // position
  position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2 * 200.0f));

  // speed and heading
  heading = glm::vec2(.0f, -1.0f);
  speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

  id = animatID;
  target = NULL;
  huntCount = 0;
  handling = false;
  handlingTimer = 0;
  wanderingTime = AppSettings::handlingTime;

  // new variables
  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  currentAttackTime = 0;
  step = 0;

  //random weights
  //ce nas zanima evolucija plena vzamemo fiksne parametre, sicer BP params
  if (PREY_EVOLUTION == 0) {

	  float sum = .0f;
	  nearestWeight = randomFloat(.0f, 1.0f);
	  sum += nearestWeight;
	  peripheralWeight = randomFloat(.0f, 1.0f);
	  sum += peripheralWeight;
	  centreWeight = randomFloat(.0f, 1.0f);
	  sum += centreWeight;

	  nearestWeight /= sum;
	  peripheralWeight /= sum;
	  centreWeight /= sum;

	  /*// random parametri z nabora moznih vrednosti (tudi v part 1)
	  distanceForAcceleration = randomFloat(1.0f, AppSettings::huntSize);
	  //velocityMultiplier = randomFloat(1.0f, 3.0f);
	  attackPeriod = randomFloat(1.0f, AppSettings::noOfSteps);*/

	  // za part 3, ko ima plenilec energijo
	  if (PRED_ENERGY_PARAMS == 1) {
		  wb = randomFloat(0.0f, 1.0f);
		  wr = randomFloat(0.0f, 1.0f);
		  we = randomFloat(0.0f, 1.0f);
		  cn0 = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);
		  cn2 = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);
		  cn1 = randomFloat(0.0f, 50.0f);
		  cn3 = randomFloat(0.0f, 50.0f);
		  cr0 = randomFloat(0.0f, energy);
		  cr1 = randomFloat(0.0f, 50.0f);

		  attackPeriod = randomInt(1, AppSettings::noOfSteps);
		  wanderingTime = randomInt(0, AppSettings::noOfSteps);
		  //velocityMultiplier = randomFloat(1.0f, 3.0f);
		  distanceForAcceleration = randomFloat(1.0f, AppSettings::huntSize);
		  //nextDecision = restPeriod;

		  v_b = randomFloat(AppSettings::cruisingSpeed, AppSettings::maxPredatorVelocity);
		  v_r = randomFloat(AppSettings::minPredatorVelocity, AppSettings::cruisingSpeed);
	  }
  }
  // v part 2 ima vse parametre fiksno dolocene, te 3 k jih ma
  else {

	  // for conf = 0 and conf = 2

	  nearestWeight = 1.0f;
	  centreWeight = .0f;
	  peripheralWeight = 0.0f;

	  // end for conf = 0
  }
  
}

Predator::Predator(int animatID, float nW, float pW, float cW, float dA, int aP, int rP, float wb1, float wr1, float we1, float cn01, float cn11, float cn21, float cn31, float cr01, float cr11, float vb, float vr) {
  // acceleration
  acceleration = glm::vec2(.0f, .0f);

  // position
  position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2 * 200.0f));

  // speed and heading
  heading = glm::vec2(.0f, -1.0f);
  speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

  id = animatID;
  target = NULL;
  huntCount = 0;
  handling = false;
  handlingTimer = 0;

  nearestWeight = nW;
  peripheralWeight = pW;
  centreWeight = cW;

  // new variables
  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  distanceForAcceleration = dA;
  //velocityMultiplier = vM;
  attackPeriod = aP;
  currentAttackTime = 0;

  wb = wb1;
  wr = wr1;
  we = we1;
  cn0 = cn01;
  cn2 = cn21;
  cn1 = cn11;
  cn3 = cn31;
  cr0 = cr01;
  cr1 = cr11;

  wanderingTime = rP;

  v_b = vb;
  v_r = vr;
}

void Predator::reset() {
	// acceleration
	acceleration = glm::vec2(.0f, .0f);

	// position
	position = glm::vec2(AppSettings::worldCentre.x, AppSettings::worldCentre.y + (AppSettings::preySize * 2 * 200.0f));

	// speed and heading
	heading = glm::vec2(.0f, -1.0f);
	speed = randomFloat(AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

	target = NULL;
	handling = false;
	wandering = false;
	handlingTimer = 0;
	wanderingTimer = 0;

	energy = 1.0f;
}

void Predator::calculate(std::vector<Prey>& preyAnimats) {

  // reset accelertion to 0 each cycle
  acceleration = glm::vec2(.0f, .0f);

  // only update stuff if not handling
  if (!(handling || wandering)) {
    
	  // neighbour data
    glm::vec2 huntVector = glm::vec2(.0f, .0f);
	currentAttackTime += 1;

	// za part 3
	#if (PRED_ENERGY_PARAMS == 1)
		if (currentAttackTime > attackPeriod) {
			target = NULL;
			wandering = true;
			currentAttackTime = 0;
		}
		if (target != NULL) {
			for (Prey p : preyAnimats) {
				if (p.isTarget && !p.isDead) {
					float distFromTarget = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
					if (distFromTarget < distanceForAcceleration) isNearCatch = true;
					else isNearCatch = false;
				}
			}
		}
	#endif

    // if has target
	if (target != NULL && !preyAnimats[target_id].isDead) {

		Prey& p = preyAnimats[target_id];

		// distance to predator's target
		float distFromPrey = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
		if (distFromPrey < .0f) distFromPrey = .0f;

		// catch attempt if close enough, always changes the target if this is true
		if (distFromPrey < AppSettings::catchDistance) {

			// confused method 1 (is close enough)
			#if (CONFUSABILITY == 1) 

				// number of confusors
				int noOfConfusors = 0;
				for (Prey p2 : preyAnimats) {
					float dist = glm::distance(p2.position, position) - AppSettings::preySize - AppSettings::predatorSize;
					if (dist < .0f) dist = .0f;
					if (dist < AppSettings::confusabilitySize) noOfConfusors++;
				}

				float confusability = 0;
				if(noOfConfusors > 0) confusability = 1 / (float)noOfConfusors;

				// not confused
				float random = randomFloat(.0f, 1.0f);
				if (random < confusability) {
					p.isDead = true;
					huntCount += 1;
					currentAttackTime = 0;
					handling = true;
				}
				// confused and starts to wander
				else {
					wandering = true;
				}
			#endif

			// confusability 0 or 2 he catches the target normally
			#if (CONFUSABILITY != 1) 
				p.isDead = true;
				huntCount += 1;
				currentAttackTime = 0;
				handling = true;
			#endif

			// it either catches it or changes targets
			p.isTarget = false;
			target = NULL;
			target_id = -1;
			isNearCatch = false;

		} 
		// end catch distance

		// if it's not close enough to catch the target but we are using Zheng's method of confusability
		// the predator switches between all the targets in his field of view
		#if (CONFUSABILITY == 2) 
			
			huntVector = glm::normalize(p.position - position);
			std::vector<int> targetsInAttackZone;

			int index = 0;
			glm::vec2 p_huntVector;
			float p_cos_alpha;
			float p_dist;

			for (Prey pr : preyAnimats) {
				if (!pr.isDead) {
					p_huntVector = glm::normalize(pr.position - position);
					p_cos_alpha = glm::dot(p_huntVector, huntVector);

					if (p_cos_alpha > cos(attackZoneAngle / 2.0f)) {
						p_dist = glm::distance(p.position, pr.position) - AppSettings::preySize * 2;
						if (p_dist < AppSettings::confusabilitySize / 2) {
							targetsInAttackZone.push_back(index);
						}
					}
				}
				index++;
			}

			int n = targetsInAttackZone.size();
			if (n > 0) {
				target_id = preyAnimats[std::rand() % (n)].id;
				target = &preyAnimats[target_id];
			}

		#endif

		// target out of range
		if (distFromPrey > AppSettings::huntSize) {
			p.isTarget = false;
			target = NULL;
			target_id = -1;
			currentAttackTime = 0;
			if(PRED_ENERGY_PARAMS == 1) wandering = true;
		}
    }

    // else find target if he doesn't have one
    else {
      currentTactic = selectTactic();
      // nearest by distance
      if (currentTactic == NEAREST) {
        float minDist = std::numeric_limits<float>::max();
        for(Prey p : preyAnimats) {
			float distFromPrey = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
			if (distFromPrey < .0f) distFromPrey = .0f;
			if (!p.isDead && distFromPrey < minDist && distFromPrey < AppSettings::huntSize) {
				minDist = distFromPrey;
				target = &(p);
				target_id = p.id;
			}
        }
      }
      // peripheral is the one with the largest peripheriality
      else if (currentTactic == PERIPHERAL) {
        float maxPer = std::numeric_limits<float>::min();
        for(Prey p : preyAnimats) {
			float distFromPrey = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
			if (distFromPrey < .0f) distFromPrey = .0f;
			if (!p.isDead && p.peripherality > maxPer && distFromPrey < AppSettings::huntSize && isOnFrontSideOfSchool(p)) {
				maxPer = p.peripherality;
				target = &(p);
				target_id = p.id;
            }
        }
      }
      // center is the one with the smallest peripheriality
      else if (currentTactic == CENTRE) {
        float minPer = std::numeric_limits<float>::max();
        for(Prey p : preyAnimats) {
			float distFromPrey = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
			if (distFromPrey < .0f) distFromPrey = .0f;
			if (!p.isDead && p.peripherality < minPer && distFromPrey < AppSettings::huntSize) {
				minPer = p.peripherality;
				target = &(p);
				target_id = p.id;
            }
        }
      }
    }

    // normalize and multiply with weight, it should have a target at this point
    if (target != NULL && target_id != -1) {

      huntVector = glm::normalize(preyAnimats[target_id].position - position);

	  if (!isExhausted) {

		  // for part 1 we want to ignore energy consumption, therefore maxForce
		  if (PRED_ENERGY_PARAMS == 0) {
			  acceleration = huntVector * AppSettings::maxPredatorForce;
		  }
		  // else we take the base huntVector and add the new drives
		  else {
			  acceleration = huntVector;
			  //  burst -> cruise
			  if (speed > cn0)
				  acceleration += pow(std::min(std::max((speed - cn0) / (AppSettings::maxPredatorVelocity - cn0), 0.0f), 1.0f), cn1) * glm::normalize(-getVelocity());
			  // rest -> crusie
			  if (speed < cn2){
				  float wn = pow(std::min(std::max((cn2 - speed) / (cn2 - AppSettings::minPredatorVelocity), 0.0f), 1.0f), cn3);
				  acceleration += wn * glm::normalize(getVelocity());
			  }
			  if (energy < cr0){
				  acceleration += we * std::min(std::max((cr0 - energy / cr0), 0.0f), cr1) * glm::normalize(-getVelocity());
			  }
		  }
	   }
    }
  }
  // handling or wandering
  else {

	    isNearCatch = false;

		handlingTimer++;
		if (handlingTimer > AppSettings::handlingTime) {
			handling = false;
			handlingTimer = 0;
			target = NULL;
			target_id = -1;
			currentAttackTime = 0;
		}

		wanderingTimer++;
		// ta je dolocen skozi evolucijo
		if (wanderingTimer > wanderingTime) {
			wandering = false;
			wanderingTimer = 0;
			target = NULL;
			target_id = -1;
			currentAttackTime = 0;
		}

		// wandering part
		nextDecision--;
		if (nextDecision == 0) {
			glm::vec2 direction = glm::vec2(randomFloat(-100.0f, 100.0f), randomFloat(-100.0f, 100.0f));
			glm::vec2 unit_temp = glm::normalize(direction);
			nextDecision = randomInt(10, 30);

			float p_cos_alpha = glm::dot(heading, unit_temp);

			// je znotraj pahljace
			if (p_cos_alpha > cos(std::_Pi / 2)) {
				unit = unit_temp;
			}
			else unit = -unit_temp;
		}

		if (!isExhausted) {

			acceleration = unit;
			//  burst -> cruise
			if (speed > cn0) acceleration += pow(std::min(std::max((speed - cn0) / (AppSettings::maxPredatorVelocity - cn0), 0.0f), 1.0f), cn1) * glm::normalize(-getVelocity());
			// rest -> crusie
			if (speed < cn2){
				float wn = pow(std::min(std::max((cn2 - speed) / (cn2 - AppSettings::minPredatorVelocity), 0.0f), 1.0f), cn3);
				acceleration += wn * glm::normalize(getVelocity());
			}
			if (energy < cr0) acceleration += we * std::min(std::max((cr0 - energy / cr0), 0.0f), cr1) * glm::normalize(-getVelocity());
		}
	}
}


void Predator::update(std::vector<Prey>& preyAnimats) {

	// update speed and heading
	glm::vec2 velocity = getVelocity();

	// speed depends on the phase 
	if (ENERGY > 0) {
		if (isNearCatch && (PRED_ENERGY_PARAMS == 1)) speed = v_b;
		else if (wandering) speed = v_r;
		else speed = AppSettings::cruisingSpeed;
	}
	
	// if fish has energy and we're using energy mode
	// changes the speed to 0 if the predator has no energy left
	if (PRED_ENERGY_PARAMS == 1){
		// added acceleration if fish still has energy
		if (energy <= 0) {
			isExhausted = true;
			speed = AppSettings::minPredatorVelocity; // minimum speed
		}
		else {
			step++;
		}
	}

	velocity = getVelocity();

	if (!isExhausted) velocity += acceleration;

    speed = .0f;
    float speed2 = glm::length2(velocity);
    if (speed2 > .0f) {
		speed = std::sqrt(speed2);
		heading = velocity / speed;
    }

    // limit speed
    speed = glm::clamp(speed, AppSettings::minPredatorVelocity, AppSettings::maxPredatorVelocity);

    // if we're using oxygen consimption
    #if (ENERGY == 1) 
	  oxygenConsumption = pow(62.9, (0.21 * (speed / AppSettings::predatorSize))) / 36; //mgO2 / kg h
	  oxygenConsumption = oxygenConsumption / 360;

	  if (speed > AppSettings::minPredatorVelocity) {
		  energy -= oxygenConsumption;
	  }
	  else {
		  energy += regenerationGain;
	  }
    #endif
  // if we're using Zheng's method
    #if (ENERGY == 2) 
	  prevAngle_t = angle_t;
	  angle_t = atan(heading.y / heading.x);
	  float average_speed = (AppSettings::minPredatorVelocity + AppSettings::maxPredatorVelocity) / 2;
	  float energyChange = 0.001f * (1 / average_speed) * (AppSettings::cruisingSpeed - speed) - 0.002f * (1 / pow(std::_Pi, 2)) * pow((angle_t - prevAngle_t), 2);

	  if (speed > AppSettings::cruisingSpeed) {
		  energy += energyChange;
	  }
	  else {
		  float speedPercentage = 1 - ((speed - AppSettings::minPredatorVelocity) / (AppSettings::cruisingSpeed - AppSettings::minPredatorVelocity));
		  energy += regenerationGain * speedPercentage;
	  }
    #endif

    // half regenerated
    if (isExhausted && (energy >= 0.5)) {
	    isExhausted = false;
    }

    if (energy > 1.0f) energy = 1.0f;
    position += getVelocity();
}

bool Predator::isOnFrontSideOfSchool(Prey const& prey)
{
  if (prey.peripherality == std::numeric_limits<float>::max()) return true;
  glm::vec2 a = glm::normalize(prey.position - position);
  return glm::dot(prey.peripheralityDir, a) > AppSettings::oppositeThreshold;
}

int Predator::selectTactic()
{

	if (PREY_EVOLUTION == 0) {
		float random = randomFloat(.0f, 1.0f);

		float threshold = nearestWeight;
		if (random < threshold)
			return NEAREST;

		threshold += peripheralWeight;
		if (random < threshold)
			return PERIPHERAL;

		return CENTRE;
	}
	else return NEAREST;
  
}