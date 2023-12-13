/***********************************************/
/**
* @file parallelCluster.cpp
*
* @brief Wrapper for Message Passing Interface (MPI).
*
* @author Torsten Mayer-Guerr
* @date 2004-11-13
*
*/
/***********************************************/

#include "base/import.h"
#include "base/string.h"
#include "base/gnssType.h"
#include "parallel/parallel.h"

// #undef SEEK_SET
// #undef SEEK_END
// #undef SEEK_CUR
#include <mpi.h>

/***** FUNCTIONS *******************************/

namespace Parallel
{

/***********************************************/

inline void check(int errorcode)
{
  if(errorcode != MPI_SUCCESS)
  {
    int  resultlen;
    char errorStr[MPI_MAX_ERROR_STRING];
    MPI_Error_string(errorcode, errorStr, &resultlen);
    throw(Exception("MPI Error: "+std::string(errorStr, resultlen)));
  }
}

/***********************************************/
/***********************************************/

// Extra channels to communicate log messages and exceptions
class HiddenChannel
{
public:
  MPI_Comm    comm;
  MPI_Request request;

  HiddenChannel(MPI_Comm communicator)
  {
    try
    {
      request = MPI_REQUEST_NULL;
      check(MPI_Comm_dup(communicator, &comm));
    }
    catch(std::exception &e)
    {
      GROOPS_RETHROW(e)
    }
  }

  virtual ~HiddenChannel()
  {
    int completed;
    check(MPI_Test(&request, &completed, MPI_STATUS_IGNORE));
    if(!completed)
    {
      check(MPI_Cancel(&request));
      check(MPI_Request_free(&request));
    }
    check(MPI_Comm_free(&comm));
  }

  virtual void receivedSignal() = 0;
};

/***********************************************/
/***********************************************/

class Mpi
{
public:
  std::vector<std::shared_ptr<HiddenChannel>> channels;

  Mpi(int argc, char *argv[])
  {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  }

 ~Mpi()
  {
    channels.clear();
    MPI_Finalize();
  }
};


/***********************************************/
/***********************************************/

class Communicator
{
public:
  std::shared_ptr<Mpi> mpi;
  MPI_Comm             comm;

  Communicator(CommunicatorPtr commParent, MPI_Comm comm_) : comm(comm_)
  {
    if(commParent)
      mpi = commParent->mpi;
  }

 ~Communicator()
  {
    if((comm != MPI_COMM_WORLD) && (comm != MPI_COMM_SELF) && (comm != MPI_COMM_NULL))
      MPI_Comm_free(&comm);
  }

  void peek();
  void wait(MPI_Request &request);
};

/***********************************************/

inline void Communicator::peek()
{
  try
  {
    if(mpi)
      for(auto &channel : mpi->channels)
        if(channel->request != MPI_REQUEST_NULL)
        {
          int flag = 0;
          check(MPI_Test(&(channel->request), &flag, MPI_STATUS_IGNORE));
          if(flag)
            channel->receivedSignal();
        }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void Communicator::wait(MPI_Request &request)
{
  try
  {
    if(request == MPI_REQUEST_NULL)
    {
      peek();
      return;
    }

    // repeat until our request finished
    for(;;)
    {
      std::vector<MPI_Request> requests(1, request);
      if(mpi)
        for(auto &channel : mpi->channels)
          requests.push_back(channel->request);

      int index = -1;
      check(MPI_Waitany(requests.size(), requests.data(), &index, MPI_STATUS_IGNORE));
      if((index == MPI_UNDEFINED) || (index < 0))
        throw(Exception("no active handles"));

      // our is request finished
      if(index == 0)
      {
        request = requests.at(0); // should be MPI_REQUEST_NULL
        peek();
        break;
      }

      // check extra channels
      mpi->channels.at(index-1)->request = requests.at(index); // should be MPI_REQUEST_NULL
      mpi->channels.at(index-1)->receivedSignal();
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

CommunicatorPtr init(int argc, char *argv[])
{
  try
  {
    auto mpi = std::make_shared<Mpi>(argc, argv);
    CommunicatorPtr comm = std::make_shared<Communicator>(nullptr, MPI_COMM_WORLD);
    comm->mpi = mpi;
    return comm;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr splitCommunicator(UInt color, UInt key, CommunicatorPtr comm)
{
  try
  {
    barrier(comm); // check for possible exceptions
    MPI_Comm commNew;
    check(MPI_Comm_split(comm->comm, ((color==NULLINDEX) ? MPI_UNDEFINED : color), key, &commNew));
    if(color==NULLINDEX)
      return nullptr;
    return std::make_shared<Communicator>(comm, commNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr createCommunicator(std::vector<UInt> ranks, CommunicatorPtr comm)
{
  try
  {
    std::vector<int> mpiRanks(ranks.size());
    for(UInt i=0; i<ranks.size(); i++)
      mpiRanks.at(i) = static_cast<int>(ranks.at(i));

    MPI_Group group, newgroup;
    check(MPI_Comm_group(comm->comm, &group));
    check(MPI_Group_incl(group, mpiRanks.size(), mpiRanks.data(), &newgroup));
#if MPI_VERSION >= 3
    // check for possible exceptions
    UInt dummy = 0;
    if(ranks.at(0) == myRank(comm))
    {
      for(UInt i=1; i<ranks.size(); i++)
        receive(dummy, ranks.at(i), comm);
      for(UInt i=1; i<ranks.size(); i++)
        send(dummy, ranks.at(i), comm);
    }
    else
      for(UInt i=1; i<ranks.size(); i++)
        if(ranks.at(i) == myRank(comm))
        {
          send   (dummy, ranks.at(0), comm);
          receive(dummy, ranks.at(0), comm);
        }
    MPI_Comm commNew;
    check(MPI_Comm_create_group(comm->comm, newgroup, 99, &commNew));
#else
    barrier(comm); // check for possible exceptions
    MPI_Comm commNew;
    check(MPI_Comm_create(comm->comm, newgroup, &commNew));
#endif
    check(MPI_Group_free(&group));
    check(MPI_Group_free(&newgroup));
    if(commNew == MPI_COMM_NULL)
      return nullptr;
    return std::make_shared<Communicator>(comm, commNew);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

CommunicatorPtr selfCommunicator()
{
  try
  {
    return std::make_shared<Communicator>(nullptr, MPI_COMM_SELF);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

UInt myRank(CommunicatorPtr comm)
{
  int rank;
  check(MPI_Comm_rank(comm->comm, &rank));
  return static_cast<UInt>(rank);
}

/***********************************************/

UInt size(CommunicatorPtr comm)
{
  int size;
  check(MPI_Comm_size(comm->comm, &size));
  return static_cast<UInt>(size);
}

/***********************************************/

void barrier(CommunicatorPtr comm)
{
  try
  {
    MPI_Request request;
    check(MPI_Ibarrier(comm->comm, &request));
    comm->wait(request);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void peek(CommunicatorPtr comm)
{
  try
  {
    comm->peek();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

class LogChannel : public HiddenChannel
{
  int                      rank, sendRank;
  unsigned int             type, count;
  std::vector<char>        data;
  std::vector<MPI_Request> sendRequests;
  std::function<void(UInt rank, UInt type, const std::string &str)> receive;

public:
  LogChannel(MPI_Comm communicator, std::function<void(UInt rank, UInt type, const std::string &str)> recv);
 ~LogChannel();
  void receivedSignal() override;
  void sendSignal(UInt type, const std::string &str);
};

/***********************************************/

LogChannel::LogChannel(MPI_Comm communicator, std::function<void(UInt rank, UInt type, const std::string &str)> recv) : HiddenChannel(communicator), sendRequests(4, MPI_REQUEST_NULL), receive(recv)
{
  try
  {
    check(MPI_Comm_rank(comm, &rank));
    if(rank == 0)
      check(MPI_Irecv(&sendRank, 1, MPI_INT, MPI_ANY_SOURCE, 111, comm, &request));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

LogChannel::~LogChannel()
{
  for(MPI_Request sendRequest : sendRequests)
  {
    int completed;
    check(MPI_Test(&sendRequest, &completed, MPI_STATUS_IGNORE));
    if(!completed)
    {
      check(MPI_Cancel(&sendRequest));
      check(MPI_Request_free(&sendRequest));
    }
  }
}

/***********************************************/

void LogChannel::sendSignal(UInt type_, const std::string &str)
{
  try
  {
    if(rank == 0)
    {
      receive(0, type_, str);
      return;
    }

    // wait for pending sends
    check(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE));

    type  = static_cast<unsigned int>(type_);
    count = str.size();
    data.resize(count);
    std::copy(str.begin(), str.end(), data.begin());

    check(MPI_Isend(&rank,  1, MPI_INT,      0, 111, comm, &sendRequests.at(0)));
    check(MPI_Isend(&type,  1, MPI_UNSIGNED, 0, 222, comm, &sendRequests.at(1)));
    check(MPI_Isend(&count, 1, MPI_UNSIGNED, 0, 333, comm, &sendRequests.at(2)));
    if(count)
      check(MPI_Isend(data.data(), count, MPI_CHAR, 0, 444, comm, &sendRequests.at(3)));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void LogChannel::receivedSignal()
{
  try
  {
    int completed = 0;
    do
    {
      check(MPI_Recv(&type,  1, MPI_UNSIGNED, sendRank, 222, comm, MPI_STATUS_IGNORE));
      check(MPI_Recv(&count, 1, MPI_UNSIGNED, sendRank, 333, comm, MPI_STATUS_IGNORE));
      data.resize(count);
      if(count)
        check(MPI_Recv(data.data(), count, MPI_CHAR, sendRank, 444, comm, MPI_STATUS_IGNORE));
      receive(static_cast<UInt>(sendRank), static_cast<UInt>(type), std::string(data.data(), count));
      data.clear();
      check(MPI_Irecv(&sendRank, 1, MPI_INT, MPI_ANY_SOURCE, 111, comm, &request));
      check(MPI_Test(&request, &completed, MPI_STATUS_IGNORE));
    }
    while(completed);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

std::function<void(UInt type, const std::string &str)> addChannel(const std::function<void(UInt rank, UInt type, const std::string &str)> &receive, CommunicatorPtr comm)
{
  try
  {
    auto channel = std::make_shared<LogChannel>(MPI_COMM_WORLD, receive);
    comm->mpi->channels.insert(comm->mpi->channels.begin(), channel); // handle logs before exceptions
    return std::bind(&LogChannel::sendSignal, channel.get(), std::placeholders::_1, std::placeholders::_2);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

class ErrorChannel : public HiddenChannel
{
  std::shared_ptr<Mpi> mpi;
  int                  code;
  unsigned int         isExternalSignal;

public:
  ErrorChannel(CommunicatorPtr comm);
  void receivedSignal() override;
  void broadCastException(const std::string &msg);
  void synchronizeAndThrowException(const std::string &msg);
};

/***********************************************/

class BroadCastedException
{
public:
  std::string message;
  BroadCastedException(const std::string &msg)  : message(msg) {}
};

static const std::string STRING_EXCEPTION_EXTERNAL("interrupted by exception from other process");

/***********************************************/

ErrorChannel::ErrorChannel(CommunicatorPtr comm)  : HiddenChannel(comm->comm), mpi(comm->mpi), code(0), isExternalSignal(0)
{
  try
  {
    check(MPI_Irecv(&code, 1, MPI_INT, MPI_ANY_SOURCE, 666, this->comm, &request));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ErrorChannel::receivedSignal()
{
  // mark all subgroups
  Bool isSubGroup = FALSE;
  auto iter = std::find_if(mpi->channels.begin(), mpi->channels.end(), [this](auto &x) {return x.get() == this;});
  if(iter != mpi->channels.end())
    while(++iter != mpi->channels.end())
    {
      if(dynamic_cast<ErrorChannel*>(iter->get()))
        dynamic_cast<ErrorChannel*>(iter->get())->isExternalSignal = 1;
      isSubGroup = TRUE;
    }

  if(isSubGroup)
    throw(Exception(STRING_EXCEPTION_EXTERNAL)); // handle subgroup first

  check(MPI_Barrier(comm));
  synchronizeAndThrowException(std::string());
}

/***********************************************/

void ErrorChannel::broadCastException(const std::string &msg)
{
  try
  {
    int rank, size;
    check(MPI_Comm_rank(comm, &rank));
    check(MPI_Comm_size(comm, &size));

    // send signal to all other processes
    int code;
    std::vector<MPI_Request> requests(size, MPI_REQUEST_NULL);
    for(int i=0; i<size; i++)
      if(i != rank)
        check(MPI_Issend(&code, 1, MPI_INT, i, 666, comm, &requests.at(i)));

    // wait for all processed (corresponding barrier() in receivedSignal())
    check(MPI_Barrier(comm));

    // cleanup all pending requests
    for(UInt i=0; i<requests.size(); i++)
      if(requests.at(i) != MPI_REQUEST_NULL)
      {
        int completed = 0;
        check(MPI_Test(&requests.at(i), &completed, MPI_STATUS_IGNORE));
        if(completed)
          continue;
        check(MPI_Cancel(&requests.at(i)));
        check(MPI_Request_free(&requests.at(i)));
      }

    synchronizeAndThrowException(msg);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void ErrorChannel::synchronizeAndThrowException(const std::string &msg)
{
  try
  {
    // complete remaining incoming messages
    int completed = 1;
    check(MPI_Test(&request, &completed, MPI_STATUS_IGNORE));
    while(completed)
    {
      check(MPI_Irecv(&code, 1, MPI_INT, MPI_ANY_SOURCE, 666, comm, &request));
      check(MPI_Test(&request, &completed, MPI_STATUS_IGNORE));
    }

    // check for any external signals
    unsigned int isExternalSignalSum = 0;
    check(MPI_Allreduce(&isExternalSignal, &isExternalSignalSum, 1, MPI_UNSIGNED, MPI_SUM, comm));
    if(isExternalSignalSum)
      throw(BroadCastedException(STRING_EXCEPTION_EXTERNAL));

    int rank, size;
    check(MPI_Comm_rank(comm, &rank));
    check(MPI_Comm_size(comm, &size));

    // collect and broadCast all messages
    std::vector<std::pair<UInt, std::string>> msgs;
    for(int process=0; process<size; process++)
    {
      unsigned int count = msg.size();
      check(MPI_Bcast(&count, 1, MPI_UNSIGNED, process, comm));
      if(!count)
        continue;
      std::vector<char> data(count);
      if(process == rank)
        std::copy(msg.begin(), msg.end(), data.begin());
      check(MPI_Bcast(data.data(), count, MPI_CHAR, process, comm));
      msgs.push_back(std::pair<UInt, std::string>(process, std::string(data.data(), count)));
    }

    // remove duplicates
    std::sort(msgs.begin(), msgs.end(), [](const auto &x, const auto &y){return (x.second == y.second) ? (x.first<y.first) : (x.second<y.second);});
    msgs.erase(std::unique(msgs.begin(), msgs.end(), [](const auto &x, const auto &y){return (x.second == y.second);}), msgs.end());
    std::sort(msgs.begin(), msgs.end(), [](const auto &x, const auto &y){return (x.first<y.first);});

    // collect messages
    std::string completeMsg;
    Bool newLine = FALSE;
    for(auto &msg : msgs)
      if(msg.second.find(STRING_EXCEPTION_EXTERNAL) == std::string::npos)
      {
        if(newLine)
          completeMsg += "\n";
        newLine = TRUE;
        if(msg.first)
          completeMsg += "in process "+msg.first%"%i of "s+size%"%i:\n"s;
        completeMsg += msg.second;
      }
    if(completeMsg.empty())
      throw(BroadCastedException("no message"));
    throw(BroadCastedException(completeMsg));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void broadCastExceptions(CommunicatorPtr commOld, std::function<void(CommunicatorPtr)> func)
{
  try // distributed exceptions
  {
    // duplicate communicator
    // as comm might be unusable after exception due to unfinished communications
    if(!commOld->mpi)
      throw(Exception("self communicator not allowed"));
    MPI_Comm commDuplicate;
    check(MPI_Comm_dup(commOld->comm, &commDuplicate));
    CommunicatorPtr comm = std::make_shared<Communicator>(commOld, commDuplicate);

    // and setup error channel
    auto errorChannel = std::make_shared<ErrorChannel>(comm);
    comm->mpi->channels.push_back(errorChannel);

    try // local exceptions
    {
      func(comm);
      barrier(comm);
    }
    catch(BroadCastedException &e)
    {
      throw;
    }
    catch(std::exception &e)
    {
      errorChannel->broadCastException(e.what());
    }
    catch(...)
    {
      errorChannel->broadCastException("Unknown ERROR");
    }

    comm->mpi->channels.pop_back();
  }
  catch(BroadCastedException &e)
  {
    commOld->mpi->channels.pop_back();
    throw(Exception(e.message));
  }
  catch(std::exception &e)
  {
    std::cerr<<"\033[1;31m"<<"****** Fatal error ******"<<std::endl;
    std::cerr<<e.what()<<"\033[0m"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    exit(EXIT_FAILURE);
  }
  catch(...)
  {
    std::cerr<<"\033[1;31m"<<"****** Fatal unknown error ******"<<"\033[0m"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    exit(EXIT_FAILURE);
  }
}

/***********************************************/

Bool isExternal(std::exception &e)
{
  return e.what() == STRING_EXCEPTION_EXTERNAL;
}

/***********************************************/
/***********************************************/

inline void send(const void *buffer, UInt count, MPI_Datatype datatype, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Request request;
    check(MPI_Isend(buffer, count, datatype, process, 17, comm->comm, &request));
    comm->wait(request);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void send(const Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(x, size, MPI_CHAR, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    send(&y, 1, MPI_UNSIGNED, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    send(&y, 1, MPI_UNSIGNED, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(static_cast<Double>(x.mjdInt()), process, comm);
    send(x.mjdMod(), process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x.type);
    send(&y, 1, MPI_UNSIGNED, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Vector3d &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(x.vector().field(), 3, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    send(x.rows(), process, comm);
    if(x.size()!=0)
      send(x.field(), x.size(), MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void send(const Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt type  = static_cast<UInt>(x.getType());
    UInt upper = x.isUpper();
    send(x.rows(),    process, comm);
    send(x.columns(), process, comm);
    send(type,        process, comm);
    send(upper,       process, comm);
    if(x.size()!=0)
      send(x.field(), x.size(), MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

inline void receive(void *buffer, UInt count, MPI_Datatype datatype, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Request request;
    check(MPI_Irecv(buffer, count, datatype, ((process!=NULLINDEX) ? process : MPI_ANY_SOURCE), 17, comm->comm, &request));
    comm->wait(request);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void receive(Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    receive(x, size, MPI_CHAR, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    receive(&y, 1, MPI_UNSIGNED, process, comm);
    x = static_cast<UInt>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    receive(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    receive(&y, 1, MPI_UNSIGNED, process, comm);
    x = static_cast<Bool>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    receive(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double mjdInt, mjdMod;
    receive(mjdInt, process, comm);
    receive(mjdMod, process, comm);
    x = Time(static_cast<Int>(mjdInt), mjdMod);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Vector3d &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double v[3];
    receive(v, 3, MPI_DOUBLE, process, comm);
    x.x() = v[0];
    x.y() = v[1];
    x.z() = v[2];
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows;
    receive(rows, process, comm);
    x = Vector(rows);
    if(x.size()!=0)
      receive(x.field(), x.size(), MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows, columns;
    UInt  type, upper;
    receive(rows,    process, comm);
    receive(columns, process, comm);
    receive(type,    process, comm);
    receive(upper,   process, comm);
    x = Matrix(rows, columns);
    x.setType(static_cast<Matrix::Type>(type), (upper) ? Matrix::UPPER : Matrix::LOWER);
    if(x.size()!=0)
      receive(x.field(), x.size(), MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void receive(GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y;
    receive(&y, 1, MPI_UNSIGNED, process, comm);
    x = GnssType(static_cast<UInt>(y));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

inline void broadCast(void *buffer, UInt count, MPI_Datatype datatype, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Request request;
    check(MPI_Ibcast(buffer, count, datatype, process, comm->comm, &request));
    comm->wait(request);
    barrier(comm); // prevents rare deadlocks when broadcast is called rapidly within loop, possibly MPI issue?
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void broadCast(Byte *x, UInt size, UInt process, CommunicatorPtr comm)
{
  try
  {
    broadCast(x, size, MPI_CHAR, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    broadCast(&y, 1, MPI_UNSIGNED, process, comm);
    x = static_cast<UInt>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    broadCast(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x);
    broadCast(&y, 1, MPI_UNSIGNED, process, comm);
    x = static_cast<Bool>(y);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Angle &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    broadCast(&x, 1, MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Time &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    Double mjdInt = x.mjdInt();
    Double mjdMod = x.mjdMod();
    broadCast(mjdInt, process, comm);
    broadCast(mjdMod, process, comm);
    x = Time(static_cast<Int>(mjdInt), mjdMod);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(GnssType &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    unsigned int y = static_cast<unsigned int>(x.type);
    broadCast(&y, 1, MPI_UNSIGNED, process, comm);
    x = GnssType(static_cast<UInt>(y));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Vector &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows = x.rows();
    broadCast(rows, process, comm);
    if(x.rows() != rows)
      x = Vector(rows);
    if(x.size()!=0)
      broadCast(x.field(), x.size(), MPI_DOUBLE, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void broadCast(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    UInt rows = x.rows();
    UInt cols = x.columns();
    UInt type  = static_cast<UInt>(x.getType());
    UInt upper = x.isUpper();

    broadCast(rows,  process, comm);
    broadCast(cols,  process, comm);
    broadCast(type,  process, comm);
    broadCast(upper, process, comm);

    if((x.rows() != rows) || (x.columns() != cols))
      x = Matrix(rows, cols);
    x.setType(static_cast<Matrix::Type>(type), (upper) ? Matrix::UPPER : Matrix::LOWER);

    constexpr UInt BLOCKSIZE = 200*1024*1024/sizeof(Double); // 200 Mb
    UInt index = 0;
    while(index<x.size())
    {
      const UInt size = std::min(x.size()-index, BLOCKSIZE);
      broadCast(x.field()+index, size, MPI_DOUBLE, process, comm);
      index += size;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

inline void reduce(const void *sendbuf, void *recvbuf, UInt count, MPI_Datatype datatype,
                   MPI_Op op, UInt process, CommunicatorPtr comm)
{
  try
  {
    MPI_Request request;
    check(MPI_Ireduce(sendbuf, recvbuf, count, datatype, op, process, comm->comm, &request));
    comm->wait(request);
    barrier(comm); // prevents rare deadlocks when reduce is called rapidly within loop, possibly MPI issue?
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int >(x);
      unsigned int y = 0;
      reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_SUM, process, comm);
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_SUM, process, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      x = 0;
      reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, process, comm);
    }
    else
      reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_SUM, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Bool &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int >(x);
      unsigned int y = 0;
      reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_SUM, process, comm);
      x = static_cast<Bool>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_SUM, process, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(Matrix &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    constexpr UInt BLOCKSIZE = 50*1024*1024/sizeof(Double); // 50 Mb

    UInt index = 0;
    while(index<x.size())
    {
      const UInt size = std::min(x.size()-index, BLOCKSIZE);
      if(myRank(comm) == process)
      {
        Vector tmp(size);
        reduce(x.field()+index, tmp.field(), size, MPI_DOUBLE, MPI_SUM, process, comm);
        std::copy_n(tmp.field(), size, x.field()+index);
      }
      else
        reduce(x.field()+index, nullptr, size, MPI_DOUBLE, MPI_SUM, process, comm);

      index += size;
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceSum(std::vector<Double> &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      std::vector<Double> tmp(x.size(), 0.);
      reduce(x.data(), tmp.data(), x.size(), MPI_DOUBLE, MPI_SUM, process, comm);
      std::swap(tmp, x);
    }
    else
      reduce(x.data(), nullptr, x.size(), MPI_DOUBLE, MPI_SUM, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void reduceMin(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int>(x);
      unsigned int y   = tmp;
      reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_MIN, process, comm);
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_MIN, process, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMin(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_MIN, process, comm);
    }
    else
      reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_MIN, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMax(UInt &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      unsigned int tmp = static_cast<unsigned int>(x);
      unsigned int y   = tmp;
      reduce(&tmp, &y, 1, MPI_UNSIGNED, MPI_MAX, process, comm);
      x = static_cast<UInt>(y);
    }
    else
    {
      unsigned int y = static_cast<unsigned int>(x);
      reduce(&y, nullptr, 1, MPI_UNSIGNED, MPI_MAX, process, comm);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void reduceMax(Double &x, UInt process, CommunicatorPtr comm)
{
  try
  {
    if(myRank(comm) == process)
    {
      Double tmp = x;
      reduce(&tmp, &x, 1, MPI_DOUBLE, MPI_MAX, process, comm);
    }
    else
      reduce(&x, nullptr, 1, MPI_DOUBLE, MPI_MAX, process, comm);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

} // namespace Parallel
