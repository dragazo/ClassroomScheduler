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

	explicit operator bool() const noexcept { return start >= 0; }
	friend bool operator<(timespan a, timespan b) noexcept { return a.start < b.start; }
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

	friend bool operator<(const timeslot &a, const timeslot &b) noexcept { return a.time < b.time; }
};
struct schedule
{
	std::string           id;          // the schedule id (name)
	std::vector<timeslot> slots;       // the time slots associated with this schedule
};
struct course_info
{
	int         capacity = 0; // maximum number of students taking the course (one section)
	std::string notes;        // notes for the course (displayed in output)

	std::set<std::string> required_attrs;     // list of all required room attributes for this course
	std::set<std::string> parallel_courses;   // list of all parallel courses (id)
	std::set<std::string> orthogonal_courses; // list of all orthogonal courses (id)
	std::set<std::string> follow_courses;     // list of all courses that must immediately follow this one
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
std::variant<std::vector<schedule>, std::string> read_schedules(const char *path)
{
	std::ifstream f{ path };
	if (!f) return std::string("failed to open ") + path + " for reading";

	std::vector<schedule> schedules;
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

		// find the schedule to put this entry into
		auto it = std::find_if(schedules.begin(), schedules.end(), [&](const schedule &i) { return i.id == str; });
		if (it == schedules.end()) // if it doesn't already exist, create it
		{
			schedules.emplace_back();
			it = schedules.end() - 1;
			it->id = std::move(str);
		}

		// insert it into the schedule in sorted order
		auto pos = std::lower_bound(it->slots.begin(), it->slots.end(), slot);
		it->slots.insert(pos, std::move(slot));
	}

	// when we're all done with that, assert some extra usage requirements
	for (const auto &schedule : schedules)
	{
		assert(!schedule.slots.empty()); // this should be impossible
		
		// make sure no timeslots overlap
		for (std::size_t i = 1; i < schedule.slots.size(); ++i)
		{
			if (schedule.slots[i].time.start < schedule.slots[i - 1].time.stop)
				return std::string(path) + " - schedule " + schedule.id + " has overlapping time slots: " + tostr(schedule.slots[i - 1].time) + " and " + tostr(schedule.slots[i].time);
		}
	}

	return std::move(schedules);
}
std::variant<std::unordered_map<std::string, course_info>, std::string> read_courses(const char *path)
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
		else return std::string(path) + ':' + std::to_string(ln) + " - unrecognized command: " + str;
	}

	return std::move(courses);
}

struct satisfy_info
{
	const std::vector<room> &rooms;
	std::vector<schedule> &schedules;
	const std::unordered_map<std::string, course_info> &courses;

	std::list<std::set<std::string>>                                            parallel_sets;          // list of parallel course sets
	std::unordered_map<std::string, std::list<std::set<std::string>>::iterator> course_to_parallel_set; // maps each course to its proper parallel set

	std::unordered_map<const std::set<std::string>*, const std::set<std::string>*> follow_psets;  // maps a pset to its follow pset (if any)
	
	std::unordered_set<const std::set<std::string>*> visited; // set of all psets that have been visited

	satisfy_info(const std::vector<room> &r, std::vector<schedule> &s, const std::unordered_map<std::string, course_info> &c) : rooms(r), schedules(s), courses(c)
	{
		// construct the trivial form of the parallel topology structure - just a bunch of bookmarked singletons
		for (const auto &entry : c)
		{
			parallel_sets.emplace_back().insert(entry.first);
			course_to_parallel_set[entry.first] = std::prev(parallel_sets.end());
		}

		// now apply all parallel constraints
		for (auto pos = parallel_sets.begin(); pos != parallel_sets.end(); )
		{
			bool merged = false;

			// for each course in the current parallel set
			for (const auto &p : *pos)
			{
				// for each parallel constraint
				for (const auto &pp : c.at(p).parallel_courses)
				{
					// get its parallel set - if it's us skip it
					auto it = course_to_parallel_set.at(pp);
					if (it == pos) continue;

					// merge it into our parallel set and account for all mapping updates
					for (const auto &ppp : *it) course_to_parallel_set[ppp] = pos;
					pos->merge(*it);
					parallel_sets.erase(it);
					merged = true;
				}
				if (merged) break; // by changing the pos set (which we're iterating over), we've invalidated the iterator and need to stop iteration (and repeat current set until no more merges)
			}

			if (!merged) ++pos; // if we didn't perform any merges we're done with this position, otherwise repeat current position
		}
		
		std::vector<std::list<std::set<std::string>>::iterator> fset;
		bool fmerge;

		// now resolve many to one follow dependencies - we can also generate the root follow chains simultaneously
		do
		{
			fmerge = false;

			for (const auto &pset : parallel_sets)
			{
				fset.clear(); // clear this for reuse (is the set of all psets that have us as a follow set)

				for (auto other = parallel_sets.begin(); other != parallel_sets.end(); ++other)
				{
					bool inserted = false;
					for (const auto &p : *other)
					{
						for (const auto &pp : c.at(p).follow_courses)
						{
							// get the follow pset
							auto val = course_to_parallel_set[pp];
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
					for (const auto &g : *fset[i]) course_to_parallel_set[g] = dest;
					dest->merge(*fset[i]);
					parallel_sets.erase(fset[i]);
				}

				// we modified the list we're iterating over and need to restart from the beginning until there are no more chances
				fmerge = true;
				break;
			}
		} while (fmerge);

		// condence heterogenous follow parallel constraints - we can also generate the follow set map simultaneously
		// this resolves one to many follow constraints
		do
		{
			fmerge = false;
			follow_psets.clear(); // clear follow pset map to undo invalidated info from previous runs

			// for each parallel set
			for (const auto &pset : parallel_sets)
			{
				// generate the fset
				fset.clear();
				for (const auto &i : pset)
					for (const auto &j : c.at(i).follow_courses)
					{
						auto val = course_to_parallel_set.at(j);
						if (std::find(fset.begin(), fset.end(), val) == fset.end()) fset.push_back(val); // add to fset but don't take duplicates
					}
				
				// if there's not enough to merge, skip to next iteration - but if there's exactly 1 item link it in the follow psets map as well
				if (fset.size() == 0) continue;
				else if (fset.size() == 1) { follow_psets[&pset] = &*fset[0]; continue; }
				
				// merge all the fsets into one parallel set
				const auto dest = fset[0];
				for (std::size_t i = 1; i < fset.size(); ++i)
				{
					for (const auto &g : *fset[i]) course_to_parallel_set[g] = dest;
					dest->merge(*fset[i]);
					parallel_sets.erase(fset[i]);
				}

				// we modified the list we're iterating over and need to restart from the beginning until there are no more chances
				fmerge = true;
				break;
			}
		} while (fmerge);

		// now we need to find the length of all follow chains
		std::unordered_map<const std::set<std::string>*, std::size_t> pset_to_follow_chain_length;
		std::vector<const std::set<std::string>*> follow_chain;

		// for each pset
		for (auto pos = parallel_sets.begin(); pos != parallel_sets.end(); ++pos)
		{
			// if this pset has already been accounted for skip it
			if (pset_to_follow_chain_length.find(&*pos) != pset_to_follow_chain_length.end()) continue;

			follow_chain.clear(); // clear follow chain for reuse

			// generate the follow chain
			for (const std::set<std::string> *p = &*pos; p; )
			{
				// if it's already in the follow chain then this is a cyclic dependency
				if (std::find(follow_chain.begin(), follow_chain.end(), p) != follow_chain.end()) throw std::logic_error("cyclic follows dependency encountered");
				follow_chain.push_back(p);

				// get the follow pset (next in the chain)
				auto it = follow_psets.find(p);
				p = it != follow_psets.end() ? it->second : nullptr;
			}

			// update follow chain lengths
			for (std::size_t i = 0; i < follow_chain.size(); ++i)
			{
				auto old_len = pset_to_follow_chain_length[follow_chain[i]];
				pset_to_follow_chain_length[follow_chain[i]] = std::max(old_len, follow_chain.size() - i); // due to arbitrary iteration order, we just need to take the max found length for each pset
			}
		}

		// we'll traverse the psets in the same order as we exit from this constructor.
		// in order to guarantee follow constraints are always satisfied we need to make sure that follow chain roots come before other items in their chain.
		// for this it is sufficient to sort by descending follow chain length - we already know that follow chains are disjoint.
		parallel_sets.sort([&](const auto &a, const auto &b) {
			auto x = pset_to_follow_chain_length.at(&a);
			auto y = pset_to_follow_chain_length.at(&b);

			return x > y || x == y && a.size() > b.size(); // we arbitrarily break ties with pset size - not required for constraint validity, but could improve performance
		});
		
		// sanity check
		for (const auto &entry : c)
		{
			auto &set = *course_to_parallel_set.at(entry.first);
			assert(set.find(entry.first) != set.end());
		}

		std::cerr << "course topology: " << parallel_sets.size() << '\n';
		for (const auto &entry : parallel_sets)
		{
			std::cerr << "[ " << std::setw(2) << pset_to_follow_chain_length.at(&entry) << " ]";
			for (const auto &s : entry) std::cerr << ' ' << std::setw(16) << s;
			std::cerr << '\n';
		}
		std::cerr << '\n';
	}
};
bool satisfy_pset_recursive_base(satisfy_info &p, std::list<std::set<std::string>>::const_iterator pset);
bool satisfy_pset_recursive_init(satisfy_info &p, schedule &sched, const std::size_t slot_i, const std::list<std::set<std::string>>::const_iterator pset_root, const std::set<std::string> *const pset);
bool satisfy_pset_recursive_interior(satisfy_info &p, schedule &sched, const std::size_t slot_i, std::list<std::set<std::string>>::const_iterator pset_root, const std::set<std::string> *const pset, std::set<std::string>::const_iterator ppos)
{
	// if we're at the end of the current pset, we're done with this pset
	if (ppos == pset->end())
	{
		// get the follow pset - if it exists then we need to recurse to it on the next slot
		if (auto it = p.follow_psets.find(pset); it != p.follow_psets.end())
		{
			// if there is no next slot, this is a failure
			if (slot_i + 1 >= sched.slots.size()) return false;

			// otherwise recurse into the follow pset
			return satisfy_pset_recursive_init(p, sched, slot_i + 1, pset_root, it->second);
		}

		// we successfully scheduled this pset and all of its (recursive) follow constraints - on to next non-follow pset from root order (get next nonvisited)
		for (++pset_root; pset_root != p.parallel_sets.end() && p.visited.find(&*pset_root) != p.visited.end(); ++pset_root) {}
		return satisfy_pset_recursive_base(p, pset_root);
	}

	const auto &course = p.courses.at(*ppos);
	timeslot   &slot = sched.slots[slot_i];

	// attempt to put it into each available room
	for (std::size_t k = 0; k < slot.assignments.size(); ++k)
	{
		// if this room is already taken, it's not viable
		if (!slot.assignments[k].empty()) continue;
		// if this room doesn't have a high enough capacity, it's not viable
		if (p.rooms[k].capacity < course.capacity) continue;

		// if this room lacks any required attributes, it's not viable
		if ([&] {
			const auto &attrs = p.rooms[k].attr;
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

		// parallel constraints are satisfied implicitly by iterating over (transitive) parallel sets

		// ----------------------------------------------------------------

		// assign current course to this schedule, timeslot, and room
		slot.assignments[k] = *ppos;

		// perform the recursive step - if we succeed, propagate up
		if (satisfy_pset_recursive_interior(p, sched, slot_i, pset_root, pset, std::next(ppos))) return true;

		// otherwise undo the change and continue searching
		slot.assignments[k].clear();
	}

	return false;
}
bool satisfy_pset_recursive_init(satisfy_info &p, schedule &sched, const std::size_t slot_i, const std::list<std::set<std::string>>::const_iterator pset_root, const std::set<std::string> *const pset)
{
	// recurse to interior with a backtracking visited marker - if the insertion fails (already present) then it was a cyclic dependency
	if (!p.visited.insert(pset).second) return false;
	if (satisfy_pset_recursive_interior(p, sched, slot_i, pset_root, pset, pset->begin())) return true;
	p.visited.erase(pset);
	return false;
}
bool satisfy_pset_recursive_base(satisfy_info &p, std::list<std::set<std::string>>::const_iterator pset)
{
	// if we're at the end of the parallel sets, we're done and have scheduled everything (yay)
	if (pset == p.parallel_sets.end()) return true;
	
	// otherwise attempt to put current pset into each schedule
	for (auto &sched : p.schedules)
	{
		// and into each timeslot in said schedule
		for (std::size_t slot_i = 0; slot_i < sched.slots.size(); ++slot_i)
		{
			// if we can satisfy this pset with this slot
			if (satisfy_pset_recursive_init(p, sched, slot_i, pset, &*pset)) return true;
		}
	}

	// otherwise we failed to satisfy the schedule (overconstrained)
	return false;
}
bool satisfy(const std::vector<room> &rooms, std::vector<schedule> &schedules, const std::unordered_map<std::string, course_info> &courses)
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

	// create info object - if it throws that means there were impossible constraints
	std::unique_ptr<satisfy_info> p;
	try { p = std::make_unique<satisfy_info>(rooms, schedules, courses); }
	catch (const std::logic_error & e) { std::cerr << "ERROR: " << e.what() << "\n\n"; return false; }

	// recurse to satisfy all course assignments
	return satisfy_pset_recursive_base(*p, p->parallel_sets.begin());
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
			if (!asgn.empty()) ostr << fix_id(asgn) << ' ' << courses.at(asgn).notes;
		}
		ostr << " \\\\\n\\hline ";
	}

	ostr << "\n\\end{tabularx}\n\\end{table}\n";
}
void print_schedules_latex(std::ostream &ostr, const std::vector<room> &rooms, const std::vector<schedule> &schedules, const std::unordered_map<std::string, course_info> &courses)
{
	ostr << "\\documentclass[8pt]{article}\n\\usepackage[margin=0.5in]{geometry}\n\\usepackage{tabularx}\n\\begin{document}\n\n";

	for (const auto &schedule : schedules)
	{
		print_schedule_latex(ostr, rooms, schedule, courses);
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
			if (!asgn.empty()) ostr << fix_id(asgn) << ' ' << courses.at(asgn).notes;
		}
		ostr << '\n';
	}
}
void print_schedules_text(std::ostream &ostr, const std::vector<room> &rooms, const std::vector<schedule> &schedules, const std::unordered_map<std::string, course_info> &courses)
{
	for (const auto &schedule : schedules)
	{
		print_schedule_text(ostr, rooms, schedule, courses);
		ostr << '\n';
	}
}

const char *help_msg = R"(Usage: schedule [OPTION]... [rooms file] [schedule file] [courses file]
Generate class schedule given rooms, timeslots, and course info.

  --help               print this help page and exit
  --latex              generate latex source instead of tab-separated text
)";

int main(int argc, const char *const argv[])
{
	std::vector<const char*> pathspec;
	bool latex = false;

	for (int i = 1; i < argc; ++i)
	{
		if (std::strcmp(argv[i], "--help") == 0)
		{
			std::cerr << help_msg;
			return 0;
		}
		else if (std::strcmp(argv[i], "--latex") == 0) latex = true;
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

	auto _courses = read_courses(pathspec[2]);
	if (_courses.index() != 0)
	{
		std::cerr << std::get<1>(_courses) << '\n';
		return 102;
	}
	auto &courses = std::get<0>(_courses);

	// ---------------------------------------------

	// attempt to satisfy all the constraints - if we fail, print error message and exit
	if (!satisfy(rooms, schedules, courses))
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