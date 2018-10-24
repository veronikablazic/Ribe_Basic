#include "AppSettings.h"
#include "Predator.h"
#include "Prey.h"

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

  //random weights
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

  // new variables
  energy = 1;
  regenerationGain = 30.0f / 100000.0f; //regeneration in 30'

  attackZoneAngle = (std::_Pi / 6);
  angle_t = 0;

  distFromTarget = 1000000;

  distanceForAcceleration = randomFloat(1.0f, AppSettings::huntSize);
  velocityMultiplier = randomFloat(1.0f, 3.0f);
  attackPeriod = randomFloat(1.0f, AppSettings::noOfSteps);
  currentAttackTime = 0;
}

Predator::Predator(int animatID, float nW, float pW, float cW, float dA, float vM, int aP) {
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
  velocityMultiplier = vM;
  attackPeriod = aP;
  currentAttackTime = 0;

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
	handlingTimer = 0;
}

void Predator::calculate(std::vector<Prey>& preyAnimats) {
  // reset accelertion to 0 each cycle
  acceleration = glm::vec2(.0f, .0f);

  // only update stuff if not handling
  if (!handling) {
    // neighbour data
    glm::vec2 huntVector = glm::vec2(.0f, .0f);
	currentAttackTime += 1;

	#if (EVOL_PARAMETERS == 1)
		if (currentAttackTime > attackPeriod) {
			target = NULL;
			handling = true;
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
    if (target != NULL) {

      
		// check if caught anything or target out of sight
		for(Prey p : preyAnimats) {
			
			// is current target
			if ((p.id == target_id) && !p.isDead) {
				float distFromPrey = glm::distance(p.position, position) - AppSettings::preySize - AppSettings::predatorSize;
				if (distFromPrey < .0f) distFromPrey = .0f;
				
				// catch attempt if close enough
				if (distFromPrey < AppSettings::catchDistance && !p.isDead) {
					
					// it can get confused and not catch it
					#if (CONFUSABILITY == 1) 
						// confusability
						float random = randomFloat(.0f, 1.0f);
						int noOfConfusors = 0;
						for (Prey p2 : preyAnimats) {
							float dist = glm::distance(p2.position, position) - AppSettings::preySize - AppSettings::predatorSize;
							if (dist < .0f) dist = .0f;
							if (dist < AppSettings::confusabilitySize) noOfConfusors++;
						}

						float confusability = 0;
						if(noOfConfusors > 0) confusability = 1 / (float)noOfConfusors;

						// not confused
						if (random < confusability) {
							p.isDead = true;
							huntCount += 1;
							currentAttackTime = 0;
						}
					#endif
					#if (CONFUSABILITY < 2) 
						p.isDead = true;
						huntCount += 1;
						currentAttackTime = 0;
					#endif

					// it either catches it or changes targets
					p.isTarget = false;
					target = NULL;
					target_id = -1;
					handling = true;
				}
				// if it's not close enough to catch the target but we are using Zheng's method of confusibility
				#if (CONFUSABILITY == 2) 
					huntVector = glm::normalize(p.position - position);
					std::vector<int> targetsInAttackZone;
					std::vector<float> targetsInAttackZoneDistance;

					int index = 0;
					float min_distance = 100000.0f;
					for (Prey pr : preyAnimats) {
						if (!pr.isDead) {
							glm::vec2 p_huntVector = glm::normalize(pr.position - position);
							// angle between current target direction and animat p
							float p_alpha = std::acos((p_huntVector.x * huntVector.x + p_huntVector.y * huntVector.y));
							if (p_alpha < 0) p_alpha += 2 * std::_Pi;
							if (p_alpha < attackZoneAngle / 2) {
								float p_dist = glm::distance(pr.position, position) - AppSettings::preySize - AppSettings::predatorSize;
								if (p_dist < min_distance) min_distance = p_dist;
								targetsInAttackZone.push_back(index);
								targetsInAttackZoneDistance.push_back(p_dist);
							}
						}
						index++;
					}

					int n = targetsInAttackZone.size();
					std::vector<int> targetsInAttackZone2;
					if (n > 0) {
						for (int i = 0; i < n; i++){
							// izberemo samo animate, ki so znotraj tega polja
							if ((min_distance + 50.0f) > targetsInAttackZoneDistance[i]) {
								targetsInAttackZone2.push_back(targetsInAttackZone[i]);
							}
						}
					}
					else {
						target = NULL;
						target_id = -1;
						currentAttackTime = 0;
					}

					n = targetsInAttackZone2.size();
					if (n > 0) {
						index = targetsInAttackZone2[std::rand() % (n)];
						target = &(preyAnimats[index]);
					}
				#endif

				// target out of range
				if (distFromPrey > AppSettings::huntSize) {
					p.isTarget = false;
					target = NULL;
					target_id = -1;
					currentAttackTime = 0;
				}

				break;
			}
		}
    }
    // find target
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
      // centre is the one with the smallest peripheriality
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

    // normalize and multiply with weight
    if (target != NULL) {
      //if (!target->isTarget) target->isTarget = true;
      huntVector = glm::normalize(target->position - position);
      acceleration = huntVector * AppSettings::maxPredatorForce;
    }
  }
  // handling
  else {
    handlingTimer++;
    if (handlingTimer > AppSettings::handlingTime) {
      handling = false;
      handlingTimer = 0;
	  target = NULL;
	  target_id = -1;
	  currentAttackTime = 0;
    }
  }
}

void Predator::update() {

  // update speed and heading
  glm::vec2 velocity = getVelocity();

  // if fish has energy and we're using energy mode
  #if (ENERGY > 0) 
	  // added acceleration if fish still has energy
	  if (!((energy > 0) && (!isExhausted))) {
		  isExhausted = true;
		  speed = AppSettings::minPredatorVelocity; // minimum speed
		  velocity = getVelocity();
	  }

	  // if not near catch swimming x times faster than slowest speed
	  else if (!isNearCatch && EVOL_PARAMETERS == 1) {
		  speed = AppSettings::minPredatorVelocity * velocityMultiplier;
		  velocity = getVelocity();
	  }
  #endif

  velocity += acceleration;

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
	  float pi = 3.14159265358979f;
	  prevAngle_t = angle_t;
	  angle_t = atan(heading.y / heading.x);
	  float average_speed = (AppSettings::minPredatorVelocity + AppSettings::maxPredatorVelocity) / 2;
	  float energyChange = 0.001f * (1 / average_speed) * (AppSettings::minPredatorVelocity - speed) - 0.002f * (1 / pow(pi, 2)) * pow((angle_t - prevAngle_t), 2);

	  if (speed > AppSettings::minPredatorVelocity) {
		  energy += energyChange;
	  }
	  else {
		  energy += regenerationGain;
	  }
  #endif

  // half regenerated
  if (isExhausted && (energy >= 0.5)) {
	  isExhausted = false;
  }

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
  float random = randomFloat(.0f, 1.0f);

  float threshold = nearestWeight;
  if (random < threshold)
    return NEAREST;

  threshold += peripheralWeight;
  if (random < threshold)
    return PERIPHERAL;

  return CENTRE;
}