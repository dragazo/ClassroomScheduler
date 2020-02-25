#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <utility>
#include <fstream>
#include <cctype>
#include <limits>
#include <variant>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <optional>
#include <list>
#include <stdexcept>
#include <memory>

void print_time(std::ostream &ostr, int time)
{
	ostr << (time / 60) << ':' << std::setfill('0') << std::setw(2) << (time % 60);
}

template<typename ...Args>
std::string tostr(Args &&...args)
{
	std::ostringstream ss;
	(ss << ... << std::forward<Args>(args));
	return ss.str();
}
std::string fix_id(std::string s)
{
	for (char &ch : s) if (ch == '_') ch = ' ';
	if (auto p = s.find('*'); p != std::string::npos) s.erase(p);
	return s;
}

struct room
{
	std::string id;                       // room id (name)
	int         capacity = 0;             // maximum number of students in this room
	std::unordered_set<std::string> attr; // list of custom attributes
};
struct timespan
{
	int start = -1;
	int stop = -1;

	friend std::ostream &operator<<(std::ostream &ostr, timespan s)
	{
		print_time(ostr, s.start);
		ostr << '-';
		print_time(ostr, s.stop);
		return ostr;
	}
};
struct timeslot
{
	timespan                 time;        // the timespan denoting this timeslot
	std::vector<std::string> assignments; // the room assignments (course id or empty string for no assignment)
};
struct schedule
{
	std::string id;              // id for this schedule (name)
	std::vector<timeslot> slots; // the time slots associated with this schedule
};
struct course_info
{
	int         capacity = 0; // maximum number of students taking the course (one section)
	std::string instructor;   // the instructor for this course
	std::string notes;        // notes for the course (displayed in output)

	std::set<std::string> required_attrs;     // list of all required room attributes for this course
	std::set<std::string> parallel_courses;   // list of all parallel courses (id)
	std::set<std::string> orthogonal_courses; // list of all orthogonal courses (id)
	std::set<std::string> follow_courses;     // list of all courses that must immediately follow this one

	std::string schedule; // if present, this is the id of the schedule that this course MUST be in
};

// skips white space but stops on new line characters (extracts it).
// returns true if stopped due to new line character.
bool skipws(std::istream &istr)
{
	for (int c; (c = istr.peek()) != EOF && std::isspace((unsigned char)c); )
	{
		istr.get();
		if (c == '\n') return true;
	}
	return false;
}

// parses a 24H time at current position in stream - returns true on success
bool parse_time(std::istream &f, int &dest)
{
	unsigned t1, t2;
	if (int c = f.peek(); c == EOF || !std::isdigit((unsigned char)c)) return false;
	if (!(f >> t1) || t1 >= 24) return false;
	if (f.get() != ':') return false;
	if (int c = f.peek(); c == EOF || !std::isdigit((unsigned char)c)) return false;
	if (!(f >> t2) || t2 >= 60) return false;

	dest = (int)(60 * t1 + t2);
	return true;
}
bool parse_time_range(std::ifstream &f, timespan &span)
{
	return parse_time(f, span.start) && f.get() == '-' && parse_time(f, span.stop);
}

auto make_line_term(std::istream &f, int &line_number)
{
	return [&] {
		if (skipws(f)) // skip white space - if we hit a new line char we're done with this line
		{
			++line_number;
			return true;
		}
		if (f.peek() == '#') // if next char starts a comment, skip past comment and on to next line
		{
			f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			++line_number;
			return true;
		}
		return false;
	};
}

std::variant<std::vector<room>, std::string> read_rooms(const char *path)
{
	std::ifstream f{ path };
	if (!f) return std::string("failed to open ") + path + " for reading";

	std::vector<room> rooms;
	int line_number = 1;
	auto line_term = make_line_term(f, line_number);

	while (true)
	{
		if (line_term()) continue; // if this is line term, line was empty (skip)
		if (!f) { assert(f.eof()); break; } // this happens at eof
		const int ln = line_number;
		room &r = rooms.emplace_back();

		f >> r.id;
		assert(f); // guaranteed to succeed from above
		if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - encountered incomplete room entry";
		if (!(f >> r.capacity)) return std::string(path) + ':' + std::to_string(ln) + " - failed to parse room capacity";

		while (!line_term() && f)
		{
			std::string s;
			f >> s;
			assert(f); // guaranteed to succeed from above
			r.attr.insert(std::move(s));
		}
	}

	if (rooms.empty()) return std::string(path) + " - no rooms were defined";

	return std::move(rooms);
}
std::variant<std::list<schedule>, std::string> read_schedules(const char *path)
{
	std::ifstream f{ path };
	if (!f) return std::string("failed to open ") + path + " for reading";

	std::list<schedule> schedules;
	std::string str;
	int line_number = 1;
	auto line_term = make_line_term(f, line_number);

	// parse all the schedule data
	while (true)
	{
		if (line_term()) continue; // if this is line term, line was empty (skip)
		if (!f) { assert(f.eof()); break; } // this happens at eof
		const int ln = line_number;
		timeslot slot;

		f >> str;
		assert(f); // guaranteed to succeed from above
		if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - encountered incomplete schedule entry";
		if (!parse_time_range(f, slot.time)) return std::string(path) + ':' + std::to_string(ln) + " - failed to parse schedule time";
		if (!line_term() && f) return std::string(path) + ':' + std::to_string(ln) + " - unexpected tokens encountered after schedule time";

		// find the schedule with the specified name
		auto it = std::find_if(schedules.begin(), schedules.end(), [&](const auto &x) { return x.id == str; });
		// if it didn't exist, create it
		if (it == schedules.end())
		{
			schedules.emplace_back().id = std::move(str);
			it = std::prev(schedules.end());
		}

		// insert it into the schedule in sorted order
		auto pos = std::lower_bound(it->slots.begin(), it->slots.end(), slot, [](const auto &a, const auto &b) { return a.time.start < b.time.start; });
		it->slots.insert(pos, std::move(slot));
	}

	// when we're all done with that, assert some extra usage requirements
	for (const auto &sched : schedules)
	{
		assert(!sched.slots.empty()); // this should be impossible
		
		// make sure no timeslots overlap
		for (std::size_t i = 1; i < sched.slots.size(); ++i)
		{
			if (sched.slots[i].time.start < sched.slots[i - 1].time.stop)
				return std::string(path) + " - schedule " + sched.id + " has overlapping time slots: " + tostr(sched.slots[i - 1].time) + " and " + tostr(sched.slots[i].time);
		}
	}

	return std::move(schedules);
}
std::variant<std::unordered_map<std::string, course_info>, std::string> read_courses(const std::list<schedule> &schedules, const char *path)
{
	std::ifstream f{ path };
	if (!f) return std::string("failed to open ") + path + " for reading";

	std::unordered_map<std::string, course_info> courses;
	std::set<std::string> constraint_set;
	std::string str;
	int line_number = 1;
	auto line_term = make_line_term(f, line_number);

	// parse all the course data
	while (true)
	{
		if (line_term()) continue; // if this is line term, line was empty (skip)
		if (!f) { assert(f.eof()); break; } // this happens at eof
		const int ln = line_number;

		f >> str;
		assert(f); // guaranteed to succeed from above
		if (str == "course")
		{
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected course id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// insert a new course into the set (make sure it doesn't already exist)
			if (courses.find(str) != courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to respecify existing course: " + str;
			course_info &info = courses[std::move(str)];

			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected course capacity";
			if (!(f >> info.capacity)) return std::string(path) + ':' + std::to_string(ln) + " - failed to parse course capacity";

			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected instructor id";
			f >> info.instructor;
			assert(f); // guaranteed to succeed from above

			if (!line_term() && f) // if we have extra stuff it's notes for the course
			{
				std::getline(f, info.notes);
				assert(f); // guaranteed to succeed from above
				++line_number; // bump this up to account for getline() consuming the newline char
				if (auto p = info.notes.find('#'); p != std::string::npos) info.notes.erase(p); // remove comment from end of string (if present)
				info.notes.erase(info.notes.find_last_not_of(" \t\n\r\v\f") + 1); // trim white space from end of string
			}
		}
		else if (str == "require")
		{
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected course id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// find the course being constrained
			auto it = courses.find(str);
			if (it == courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str;

			// gather all the tokens
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected one or more room attributes";
			do
			{
				f >> str;
				assert(f); // guaranteed to succeed from above
				it->second.required_attrs.insert(std::move(str)); // add requirement (duplicates are no-op)
			} while (!line_term() && f);
		}
		else if (str == "parallel" || str == "orthogonal")
		{
			// parse the course list
			constraint_set.clear();
			while (!line_term() && f)
			{
				std::string s;
				f >> s;
				assert(f); // guaranteed to succeed from above
				if (courses.find(s) == courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + s;
				constraint_set.insert(std::move(s));
			}
			if (constraint_set.size() < 2) return std::string(path) + ':' + std::to_string(ln) + " - expected two or more course ids";

			// we handle parallel and orthogonal with same code - use this to separate the semantics
			const auto bucket = str == "parallel" ? &course_info::parallel_courses : &course_info::orthogonal_courses;

			// apply constraints
			for (auto &dest_id : constraint_set)
			{
				auto &dest = courses.at(dest_id).*bucket; // get the destination set
				for (auto &constraint : constraint_set) if (&dest_id != &constraint)
				{
					dest.insert(constraint); // add all other constraints to dest
				}
			}
		}
		else if (str == "follows")
		{
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected course id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// find the course being constrained (first)
			auto first = courses.find(str);
			if (first == courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str;

			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected a second course id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// find the course being constrained (second)
			auto second = courses.find(str);
			if (second == courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str;

			if (!line_term() && f) return std::string(path) + ':' + std::to_string(ln) + " - expected only two arguments";

			// apply the constraint
			first->second.follow_courses.insert(second->first);
		}
		else if (str == "schedule")
		{
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected course id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// find the course being constrained
			auto it = courses.find(str);
			if (it == courses.end()) return std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str;

			// get the schedule id
			if (line_term() || !f) return std::string(path) + ':' + std::to_string(ln) + " - expected a schedule id";
			f >> str;
			assert(f); // guaranteed to succeed from above

			// if there's already a schedule constraint for this course and they differ, it's a problem
			std::string &sch = it->second.schedule;
			if (!sch.empty() && sch != str) return std::string(path) + ':' + std::to_string(ln) + " - attempt to override pre-existing schedule constraint: " + sch + " -> " + str;

			// if this is an unknown course, it's a problem
			if (std::find_if(schedules.begin(), schedules.end(), [&](const auto &x) { return x.id == str; }) == schedules.end()) return std::string(path) + ':' + std::to_string(ln) + " - reference to undefined schedule: " + str;
			sch = str; // apply the constraint
		}
		else return std::string(path) + ':' + std::to_string(ln) + " - unrecognized command: " + str;
	}

	return std::move(courses);
}

struct satisfy_info
{
public:
	const std::vector<room> &rooms;
	std::list<schedule> &schedules;
	const std::unordered_map<std::string, course_info> &courses;

	std::list<std::list<std::set<std::string>>>                            follow_chains;            // list of all follow chains in the course topology
	std::unordered_map<const std::list<std::set<std::string>>*, schedule*> follow_chain_to_schedule; // maps each follow chain to its required schedule or nullptr if no specific schedule is required

	// constructs a new satisfiability info object
	// if this throws it means that the system is impossibly-constrained (as opposed to just being unsatisfiable)
	satisfy_info(const std::vector<room> &r, std::list<schedule> &s, const std::unordered_map<std::string, course_info> &c) : rooms(r), schedules(s), courses(c)
	{
		std::list<std::set<std::string>>                                               psets;          // list of all parallel course sets
		std::unordered_map<std::string, std::list<std::set<std::string>>::iterator>    course_to_pset; // maps each course to its pset
		std::unordered_map<const std::set<std::string>*, const std::set<std::string>*> follow_psets;   // maps a pset to its follow pset (if any)
		std::vector<std::list<std::set<std::string>>::iterator>                        fset;           // temporary for building follow sets
		bool                                                                           fmerge;         // flag for looping construction logic

		// construct the trivial form of the parallel topology structure - just a bunch of bookmarked singletons
		for (const auto &entry : c)
		{
			psets.emplace_back().insert(entry.first);
			course_to_pset[entry.first] = std::prev(psets.end());
		}

		// now apply all parallel constraints
		for (auto pos = psets.begin(); pos != psets.end(); )
		{
			bool merged = false;

			// for each course in the current parallel set
			for (const auto &p : *pos)
			{
				// for each parallel constraint
				for (const auto &pp : c.at(p).parallel_courses)
				{
					// get its parallel set - if it's us skip it
					auto it = course_to_pset.at(pp);
					if (it == pos) continue;

					// merge it into our parallel set and account for all mapping updates
					for (const auto &ppp : *it) course_to_pset[ppp] = pos;
					pos->merge(*it);
					psets.erase(it);
					merged = true;
				}
				if (merged) break; // by changing the pos set (which we're iterating over), we've invalidated the iterator and need to stop iteration (and repeat current set until no more merges)
			}

			if (!merged) ++pos; // if we didn't perform any merges we're done with this position, otherwise repeat current position
		}
		
		// now resolve many to one follow dependencies
		do
		{
			fmerge = false;

			for (const auto &pset : psets)
			{
				fset.clear(); // clear this for reuse (is the set of all psets that have us as a follow set)

				for (auto other = psets.begin(); other != psets.end(); ++other)
				{
					bool inserted = false;
					for (const auto &p : *other)
					{
						for (const auto &pp : c.at(p).follow_courses)
						{
							// get the follow pset
							auto val = course_to_pset[pp];
							// if it's us, add it to the fset
							if (&*val == &pset)
							{
								fset.push_back(other); // insert into fset and move on to next other pset (this ensures that fset has no duplicates)
								inserted = true;
								break;
							}
						}
						if (inserted) break; // propagate break request
					}
				}

				// if there aren't enough to merge, move on to next pset
				if (fset.size() < 2) continue;

				// merge all the fsets into one parallel set
				const auto dest = fset[0];
				for (std::size_t i = 1; i < fset.size(); ++i)
				{
					for (const auto &g : *fset[i]) course_to_pset[g] = dest;
					dest->merge(*fset[i]);
					psets.erase(fset[i]);
				}

				// we modified the list at unknown locations and need to restart from the beginning until there are no more chances
				fmerge = true;
				break;
			}
		} while (fmerge);

		// resolves one to many follow constraints - we can also generate the follow set map simultaneously
		do
		{
			fmerge = false;
			follow_psets.clear(); // clear follow pset map to undo invalidated info from previous runs

			// for each parallel set
			for (const auto &pset : psets)
			{
				// generate the fset
				fset.clear();
				for (const auto &i : pset)
				{
					for (const auto &j : c.at(i).follow_courses)
					{
						auto val = course_to_pset.at(j);
						if (std::find(fset.begin(), fset.end(), val) == fset.end()) fset.push_back(val); // add to fset but don't take duplicates
					}
				}
				
				// if there's not enough to merge, skip to next iteration - but if there's exactly 1 item link it in the follow psets map as well
				if (fset.size() == 0) continue;
				else if (fset.size() == 1) { follow_psets[&pset] = &*fset[0]; continue; }
				
				// merge all the fsets into one parallel set
				const auto dest = fset[0];
				for (std::size_t i = 1; i < fset.size(); ++i)
				{
					for (const auto &g : *fset[i]) course_to_pset[g] = dest;
					dest->merge(*fset[i]);
					psets.erase(fset[i]);
				}

				// we modified the list we're iterating over and need to restart from the beginning until there are no more chances
				fmerge = true;
				break;
			}
		} while (fmerge);

		// at this point all psets have 0 or 1 children and 0 or 1 parents - so everything is either a chain (good) or a cycle (impossible constraints)
		// we now need to generate the follow chains - start by creating singletons for all the psets that don't have any children (terminal)
		for (auto pos = psets.begin(); pos != psets.end(); ) 
		{
			// if this pset doesn't have a follow pset (leaf node)
			if (follow_psets.find(&*pos) == follow_psets.end())
			{
				// create a new follow chain and splice this node into it
				auto next = std::next(pos);
				auto &chain = follow_chains.emplace_back();
				chain.splice(chain.begin(), psets, pos);
				pos = next;
			}
			else ++pos;
		}
		// now repeatedly append to follow chains until there are no more changes
		do
		{
			fmerge = false;

			for (auto pos = psets.begin(); pos != psets.end(); )
			{
				// get the follow set for this pset - we know one exists because we already removed all the ones that didn't with the above loop
				auto follow = follow_psets.at(&*pos);
				bool inserted = false;

				// look for a follow chain to put it in (must go into the front to ensure we don't lazily permit cycles
				for (auto &chain : follow_chains) if (&chain.front() == follow)
				{
					auto next = std::next(pos);
					chain.splice(chain.begin(), psets, pos);
					pos = next;
					inserted = true; // this flag denotes we inserted the current item (just used for correct increment logic)
					fmerge = true;   // this flag denotes that we made a change at all (for repeat logic)
					break;
				}

				if (!inserted) ++pos;
			}
		} while (fmerge);
		// we repeated that logic until there were no more changes - if psets is not empty at this point then there were cyclic dependencies
		if (!psets.empty()) throw std::logic_error("cyclic follows dependency encountered");

		// for performance reasons when generating schedules it's better to traverse larger follow chains first
		follow_chains.sort([](const auto &a, const auto &b) { return a.size() > b.size(); });

		// now we need to generate the follow chain schedule requirements map
		for (const auto &chain : follow_chains)
		{
			schedule          *required_schedule = nullptr; // keep track of the required schedule for this follow chain
			const std::string *source = nullptr;            // also keep track of the source for this constraint for error messages
			
			for (const auto &pset : chain)
			{
				for (const auto &p : pset)
				{
					// get this course - if it doesn't require a specific schedule skip it
					const course_info &course = c.at(p);
					if (course.schedule.empty()) continue;

					// get the required schedule for this course
					auto sched = std::find_if(schedules.begin(), schedules.end(), [&](const auto &x) { return x.id == course.schedule; });
					assert(sched != schedules.end()); // this should have been guaranteed prior to invoking the constructor

					// if there's a required schedule for this chain we need to match it
					if (required_schedule)
					{
						if (required_schedule != &*sched) throw std::logic_error(std::string("conflicting schedule specifiers in same follow chain: ") + *source + " -> " + required_schedule->id + " vs " + p + " -> " + course.schedule);
					}
					// otherwise just mark this as the new required schedule for the chain
					else
					{
						required_schedule = &*sched;
						source = &p;
					}
				}
			}

			follow_chain_to_schedule[&chain] = required_schedule; // add this chain to the map
		}
	}

	void print_topology(std::ostream &ostr)
	{
		ostr << "course topology: " << follow_chains.size() << '\n';
		for (const auto &chain : follow_chains)
		{
			ostr << "[ " << chain.size() << " ]:\n";
			for (const auto &s : chain)
			{
				ostr << "   ";
				for (const auto &ss : s) ostr << ' ' << std::setw(16) << ss;
				ostr << '\n';
			}
		}
		ostr << '\n';
	}

	// begins by erasing the current schedules, then recursively applies constraints in order to generate a valid schedule.
	// returns true on success (and leaves the result in schedule).
	bool satisfy();

private:
	bool satisfy_fchain(std::list<std::list<std::set<std::string>>>::const_iterator chain);
	bool satisfy_pset_at(std::list<std::list<std::set<std::string>>>::const_iterator chain, schedule &sched, std::size_t slot_i, std::list<std::set<std::string>>::const_iterator pset);
	bool satisfy_pset_at_interior(std::list<std::list<std::set<std::string>>>::const_iterator chain, schedule &sched, std::size_t slot_i, std::list<std::set<std::string>>::const_iterator pset, std::set<std::string>::const_iterator ppos);
};

bool satisfy_info::satisfy_pset_at_interior(std::list<std::list<std::set<std::string>>>::const_iterator chain, schedule &sched, std::size_t slot_i, std::list<std::set<std::string>>::const_iterator pset, std::set<std::string>::const_iterator ppos)
{
	// if we're at the end of the current pset, we're done with this pset
	if (ppos == pset->end())
	{
		assert(slot_i + 1 < sched.slots.size()); // this should be guaranteed at the follow chain level

		// recurse to the next item in the follow chain
		return satisfy_pset_at(chain, sched, slot_i + 1, std::next(pset));
	}

	const auto &course = courses.at(*ppos);
	timeslot   &slot = sched.slots[slot_i];

	// attempt to put it into each available room
	for (std::size_t k = 0; k < slot.assignments.size(); ++k)
	{
		// if this room is already taken, it's not viable
		if (!slot.assignments[k].empty()) continue;
		// if this room doesn't have a high enough capacity, it's not viable
		if (rooms[k].capacity < course.capacity) continue;

		// if this room lacks any required attributes, it's not viable
		if ([&] {
			const auto &attrs = rooms[k].attr;
				for (const auto &req : course.required_attrs)
				{
					if (attrs.find(req) == attrs.end()) return true;
				}
			return false;
		}()) continue;

		// if we violate an orthogonality constraint, it's not viable
		if ([&] {
			const auto &ortho = course.orthogonal_courses;
				for (const auto &other_assignment : slot.assignments)
				{
					if (ortho.find(other_assignment) != ortho.end()) return true;
				}
			return false;
		}()) continue;

		// ----------------------------------------------------------------

		// assign current course to this schedule, timeslot, and room
		slot.assignments[k] = *ppos;

		// perform the recursive step - if we succeed, propagate up
		if (satisfy_pset_at_interior(chain, sched, slot_i, pset, std::next(ppos))) return true;

		// otherwise undo the change and continue searching
		slot.assignments[k].clear();
	}

	return false;
}
bool satisfy_info::satisfy_pset_at(std::list<std::list<std::set<std::string>>>::const_iterator chain, schedule &sched, std::size_t slot_i, std::list<std::set<std::string>>::const_iterator pset)
{
	// if we're at the end of the follow chain, we're done with this chain - recurse to the next
	if (pset == chain->end()) return satisfy_fchain(std::next(chain));
	// otherwise recurse into the pset
	else return satisfy_pset_at_interior(chain, sched, slot_i, pset, pset->begin());
}
bool satisfy_info::satisfy_fchain(std::list<std::list<std::set<std::string>>>::const_iterator chain)
{
	// if we're at the end of the follow chains we're done and have scheduled everything successfully (yay)
	if (chain == follow_chains.end()) return true;

	// create a function to attempt to place the current course into the given schedule
	auto scheduler = [&] (schedule &sched) {
		// try to put it into each timeslot in the given schedule
		for (std::size_t slot_i = 0; slot_i < sched.slots.size(); ++slot_i)
		{
			// if there's not enough room to fit the entire chain we can prune this entire execution branch
			if (slot_i + chain->size() >= sched.slots.size()) break;

			// if we can satisfy this pset with this slot
			if (satisfy_pset_at(chain, sched, slot_i, chain->begin())) return true;
		}
		return false;
	};

	// get the required schedule for this chain - if there is one, use that
	if (schedule *req = follow_chain_to_schedule.at(&*chain); req)
	{
		return scheduler(*req);
	}
	// otherwise attempt to put current follow chain into each schedule in order
	else
	{
		for (auto &sched : schedules)
		{
			if (scheduler(sched)) return true;
		}

		// otherwise we failed to satisfy the scuedule (overconstrained)
		return false;
	}
}
bool satisfy_info::satisfy()
{
	// initialize all schedule assignments to empty
	for (auto &sched : schedules)
	{
		for (auto &slot : sched.slots)
		{
			slot.assignments.clear();
			slot.assignments.resize(rooms.size());
		}
	}

	// recurse to satisfy all course assignments
	return satisfy_fchain(follow_chains.begin());
}

void print_schedule_latex(std::ostream &ostr, const std::vector<room> &rooms, const schedule &sched, const std::unordered_map<std::string, course_info> &courses)
{
	ostr << "\\begin{table}[ht!]\n\\centering\n\\begin{tabularx}{\\textwidth}{|X";
	for (std::size_t i = 0; i < rooms.size(); ++i) ostr << "|X";
	ostr << "|}\n\\hline ";

	ostr << fix_id(sched.id);
	for (const auto &room : rooms) ostr << " & " << fix_id(room.id);
	ostr << " \\\\\n\\hline ";

	for (const auto &timeslot : sched.slots)
	{
		ostr << timeslot.time;
		for (const auto &asgn : timeslot.assignments)
		{
			ostr << " & ";
			if (!asgn.empty())
			{
				const auto &c = courses.at(asgn);
				ostr << fix_id(asgn) << ' ' << fix_id(c.instructor);
				if (!c.notes.empty()) ostr << ' ' << c.notes;
			}
		}
		ostr << " \\\\\n\\hline ";
	}

	ostr << "\n\\end{tabularx}\n\\end{table}\n";
}
void print_schedules_latex(std::ostream &ostr, const std::vector<room> &rooms, const std::list<schedule> &schedules, const std::unordered_map<std::string, course_info> &courses)
{
	ostr << "\\documentclass[8pt]{article}\n\\usepackage[margin=0.5in]{geometry}\n\\usepackage{tabularx}\n\\begin{document}\n\n";

	for (const auto &sched : schedules)
	{
		print_schedule_latex(ostr, rooms, sched, courses);
		ostr << '\n';
	}

	ostr << "\\end{document}\n";
}

void print_schedule_text(std::ostream &ostr, const std::vector<room> &rooms, const schedule &sched, const std::unordered_map<std::string, course_info> &courses)
{
	ostr << fix_id(sched.id);
	for (const auto &room : rooms) ostr << '\t' << fix_id(room.id);
	ostr << '\n';

	for (const auto &timeslot : sched.slots)
	{
		ostr << timeslot.time;
		for (const auto &asgn : timeslot.assignments)
		{
			ostr << '\t';
			if (!asgn.empty())
			{
				const auto &c = courses.at(asgn);
				ostr << fix_id(asgn) << ' ' << fix_id(c.instructor);
				if (!c.notes.empty()) ostr << ' ' << c.notes;
			}
		}
		ostr << '\n';
	}
}
void print_schedules_text(std::ostream &ostr, const std::vector<room> &rooms, const std::list<schedule> &schedules, const std::unordered_map<std::string, course_info> &courses)
{
	for (const auto &sched : schedules)
	{
		print_schedule_text(ostr, rooms, sched, courses);
		ostr << '\n';
	}
}

const char *const help_msg = R"(Usage: schedule [OPTION]... [rooms file] [schedule file] [courses file]
Generate class schedule given rooms, timeslots, and course info.

  --help               print this help page and exit
  --latex              generate latex source instead of tab-separated text
  --topology           print the course topology to stderr in addition to all other work
)";

int main(int argc, const char *const argv[])
{
	std::vector<const char*> pathspec;
	bool latex = false;
	bool topology = false;

	for (int i = 1; i < argc; ++i)
	{
		if (std::strcmp(argv[i], "--help") == 0)
		{
			std::cerr << help_msg;
			return 0;
		}
		else if (std::strcmp(argv[i], "--latex") == 0) latex = true;
		else if (std::strcmp(argv[i], "--topology") == 0) topology = true;
		else pathspec.push_back(argv[i]);
	}

	if (pathspec.size() != 3)
	{
		std::cerr << "incorrect usage: see --help for info\n";
		return 1;
	}
	std::cout.fill('0');

	// ---------------------------------------------

	auto _rooms = read_rooms(pathspec[0]);
	if (_rooms.index() != 0)
	{
		std::cerr << std::get<1>(_rooms) << '\n';
		return 100;
	}
	auto &rooms = std::get<0>(_rooms);

	auto _schedules = read_schedules(pathspec[1]);
	if (_schedules.index() != 0)
	{
		std::cerr << std::get<1>(_schedules) << '\n';
		return 101;
	}
	auto &schedules = std::get<0>(_schedules);

	auto _courses = read_courses(schedules, pathspec[2]);
	if (_courses.index() != 0)
	{
		std::cerr << std::get<1>(_courses) << '\n';
		return 102;
	}
	auto &courses = std::get<0>(_courses);

	// ---------------------------------------------

	// attempt to create the solver object - if this fails it means that there were impossible constraints
	std::unique_ptr<satisfy_info> solver;
	try { solver = std::make_unique<satisfy_info>(rooms, schedules, courses); }
	catch (const std::exception & e)
	{
		std::cerr << "impossible constraints encountered: " << e.what() << '\n';
		return 240;
	}

	// if topology flag was set, print topology info to stderr
	if (topology) solver->print_topology(std::cerr);

	// attempt to satisfy all the constraints - if we fail, print error message and exit
	if (!solver->satisfy())
	{
		std::cerr << "failed to generate schedule (overconstrained)\n";
		return 200;
	}

	// ---------------------------------------------

	// print the resulting (satisfied) schedule
	auto printer = latex ? print_schedules_latex : print_schedules_text;
	printer(std::cout, rooms, schedules, courses);

	return 0;
}