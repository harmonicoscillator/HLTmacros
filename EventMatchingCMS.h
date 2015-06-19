// General Purpose event matching class.
// Takes event, run, lumi, and makes a map
//
// Instructions:
// Generally you will have two files that you need to match events for.
// First, pick one file. Loop over the file and add each event
// using the addEvent method. The entry number should be the TTree entry number.
// Once the map is completed, begin your analysis loop and retrieve
// matching entries from the first file by calling retrieveKey with
// event info from the second file.
//
// The matcher will return -1 in the event of no matching event.
// The matcher does not allow duplicate entries. Keys are removed
// from the map after retrieval, preventing duplicates on both
// sides of the transaction.

#include <map>

#ifndef EVENTMATCHINGCMS
#define EVENTMATCHINGCMS

class EventMatchingCMS{
public:
  EventMatchingCMS(){};
  ~EventMatchingCMS(){};

  bool addEvent(long long event, long long lumi, long long run, long long entryNumber);
  long long retrieveEvent(long long event, long long lumi, long long run);

private:
  long long makeKey(long long event, long long lumi, long long run);
  bool addKey(long long key, long long entryNumber);

  std::map<long long, long long> kmap;

};
#endif


long long EventMatchingCMS::makeKey(long long event, long long lumi, long long run)
{
  return (10000000000000000*run + 10000000000*lumi + event);
}

bool EventMatchingCMS::addKey(long long key, long long entryNumber)
{
  std::pair<long long, long long> p(key, entryNumber);
  return kmap.insert(p).second;
}

bool EventMatchingCMS::addEvent(long long event, long long lumi, long long run, long long entryNumber)
{
  long long key = makeKey(event, lumi, run);
  return addKey(key, entryNumber);
}

long long EventMatchingCMS::retrieveEvent(long long event, long long lumi, long long run)
{
  long long key = makeKey(event, lumi, run);

  std::map<long long, long long>::const_iterator got = kmap.find(key);
  if(got == kmap.end() ) {
    return -1;
  } else {
    kmap.erase(got);
    return got->second;
  }
}
