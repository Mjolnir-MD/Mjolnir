#ifndef MJOLNIR_SIMULATOR_BASE
#define MJOLNIR_SIMULATOR_BASE

namespace mjolnir
{

class SimulatorBase
{
  public:
    virtual ~SimulatorBase() = default;

    virtual void initialize() = 0;
    virtual bool step()       = 0;
    virtual void run()        = 0;
    virtual void finalize()   = 0;
};

}// mjolnir
#endif// MJOLNIR_SIMULATOR_BASE
