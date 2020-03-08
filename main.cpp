#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <string>
#include <cstring>
#include <deque>
#include <list>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <cctype>
#include <limits>
#include <variant>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cassert>

void print_time(std::ostream &ostr, int time)
{
	auto oldfill = ostr.fill('0');
	ostr << (time / 60) << ':' << std::setw(2) << (time % 60);
	ostr.fill(oldfill);
}

template<typename ...Args>
std::string tostr(Args &&...args)
{
	std::ostringstream ss;
	(ss << ... << std::forward<Args>(args));
	return ss.str();
}
std::string to_lower(std::string s)
{
	for (char &ch : s) ch = std::tolower((unsigned char)ch);
	return s;
}
std::string fix_id(std::string s)
{
	for (char &ch : s) if (ch == '_') ch = ' ';
	if (auto p = s.find('*'); p != std::string::npos) s.erase(p);
	return s;
}

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
struct time_union
{
	std::vector<timespan> content;

	template<bool closed>
	auto intersect_range(timespan span)
	{
		if constexpr (closed)
		{
			// equal_range() gives closed intervals with less than predicates by definition
			auto pos = std::lower_bound(content.begin(), content.end(), span.start, [](const timespan &a, int b) { return a.stop < b; });
			auto other = std::upper_bound(content.begin(), content.end(), span.stop, [](int a, const timespan &b) { return a < b.start; });
			return std::make_pair(pos, other);
		}
		else
		{
			// but by giving a less than or equal predicate we get an open interval by omitting the equal items
			auto pos = std::lower_bound(content.begin(), content.end(), span.start, [](const timespan &a, int b) { return a.stop <= b; });
			auto other = std::upper_bound(content.begin(), content.end(), span.stop, [](int a, const timespan &b) { return a <= b.start; });
			return std::make_pair(pos, other);
		}
	}

	void add(timespan span)
	{
		static_assert(std::is_same_v<std::iterator_traits<decltype(content)::iterator>::iterator_category, std::random_access_iterator_tag>);
		assert(span.start >= 0 && span.start <= span.stop);

		// construct [pos, other) - the range that this span intersects (as closed intervals)
		auto range = intersect_range<true>(span);

		// if there was no intersection, just insert into the list at appropriate location
		if (range.first == range.second) content.insert(range.first, span);
		// otherwise there was an intersection and we need to do merging logic
		else
		{
			// merge the times into the first intersection item
			range.first->start = std::min(range.first->start, span.start);
			range.first->stop = std::max(std::prev(range.second)->stop, span.stop);
			// and erase all the others (now covered by first)
			content.erase(std::next(range.first), range.second);
		}
	}
	void add(const time_union u)
	{
		for (const auto &i : u.content) add(i);
	}

	bool intersects(timespan span)
	{
		auto range = intersect_range<false>(span); // just checking for intersection uses open interval logic
		return range.first != range.second;
	}

	bool empty() const noexcept { return content.empty(); }
};

struct course_t;
struct room_t
{
	std::string id;                       // room id (name)
	int         capacity = 0;             // maximum number of students in this room
	std::unordered_set<std::string> attr; // list of custom attributes
};
struct timeslot_t
{
	timespan                   time;        // the timespan denoting this timeslot
	std::vector<const course_t*> assignments; // the room assignments (course id or empty string for no assignment)
};
struct schedule_t
{
	std::string id;              // id for this schedule (name)
	std::vector<timeslot_t> slots; // the time slots associated with this schedule
};
struct instructor_t
{
	std::string id;             // id (name) of this instructor
	time_union  unavailability; // all times when this instructor is unavailable for some reason (e.g. early morning due to long commute or other responsibilities)
};
struct course_t
{
	std::string id;                   // id (name) of this course
	int           capacity = 0;         // maximum number of students taking the course (one section)
	instructor_t *instructor = nullptr; // the instructor for this course
	std::string   notes;                // notes for the course (displayed in output)

	std::set<std::string> required_attrs; // list of all required room attributes for this course

	std::set<course_t*> parallel_courses;   // list of all parallel courses (id)
	std::set<course_t*> orthogonal_courses; // list of all orthogonal courses (id)
	std::set<course_t*> follow_courses;     // list of all courses that must immediately follow this one

	schedule_t *required_schedule = nullptr; // if present, this is the id of the schedule that this course MUST be in
};
struct constraint_set
{
	std::deque<room_t>       rooms;       // a list of all the rooms
	std::deque<schedule_t>   schedules;   // a list of all the schedules
	std::deque<instructor_t> instructors; // a list of all instructors
	std::deque<course_t>     courses;     // a list of all courses

	std::unordered_map<std::string, room_t*>       rooms_map;       // maps from room id to the room
	std::unordered_map<std::string, schedule_t*>   schedules_map;   // maps from schedule id to the schedule
	std::unordered_map<std::string, instructor_t*> instructors_map; // a list of all instructors - maps id to instructor info
	std::unordered_map<std::string, course_t*>     courses_map;     // a list of all courses - maps id to course info
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
bool parse_time_range(std::istream &f, timespan &span)
{
	return parse_time(f, span.start) && f.get() == '-' && parse_time(f, span.stop);
}

struct line_term_t
{
	std::istream &file;
	int line_number;

	bool operator()()
	{
		if (skipws(file)) // skip white space - if we hit a new line char we're done with this line
		{
			++line_number;
			return true;
		}
		if (file.peek() == '#') // if next char starts a comment, skip past comment and on to next line
		{
			file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			++line_number;
			return true;
		}
		return false;
	}
};
struct constraint_parser_pack
{
	constraint_set c;
	std::set<course_t*> course_set;
	std::string str, err;
	std::istream &f;
	const char *const path;
	int ln;
	line_term_t line_term;

	constraint_parser_pack(std::istream &file, const char *p) : f(file), path(p), ln(1), line_term{ f, 1 }  {}

	bool _parse_room()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected room id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// make sure a room with the same name doesn't already exist
		if (c.courses_map.find(str) != c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - a room with this id already exists"; return false; }

		// create the new room and link to to the lookup table
		room_t &r = c.rooms.emplace_back();
		r.id = str;
		c.rooms_map.emplace(std::move(str), &r);

		// parse capacity
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected room capacity"; return false; }
		if (!(f >> r.capacity)) { err = std::string(path) + ':' + std::to_string(ln) + " - failed to parse room capacity"; return false; }
		if (int ch = f.peek(); ch != EOF && !std::isspace((unsigned char)ch)) { err = std::string(path) + ':' + std::to_string(ln) + " - unexpected character encountered in room capacity"; return false; }

		// parse (optional) attributes
		while (!line_term() && f)
		{
			f >> str;
			assert(f); // guaranteed to succeed from above
			r.attr.insert(std::move(str));
		}

		return true;
	}
	bool _parse_timeslot()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected schedule id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// get the schedule being modified - if it doesn't already exist create it
		auto it = c.schedules_map.find(str);
		if (it == c.schedules_map.end())
		{
			schedule_t &s = c.schedules.emplace_back();
			s.id = str;
			it = c.schedules_map.emplace(std::move(str), &s).first;
		}

		timeslot_t slot;

		// parse time points
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - encountered incomplete schedule entry"; return false; }
		if (!parse_time_range(f, slot.time)) { err = std::string(path) + ':' + std::to_string(ln) + " - failed to parse schedule time"; return false; }
		if (!line_term() && f) { err = std::string(path) + ':' + std::to_string(ln) + " - unexpected tokens encountered after schedule time"; return false; }

		// insert it into the schedule in sorted order
		auto pos = std::lower_bound(it->second->slots.begin(), it->second->slots.end(), slot, [](const auto &a, const auto &b) { return a.time.start < b.time.start; });
		it->second->slots.insert(pos, std::move(slot));

		return true;
	}
	bool _parse_instructor()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected instructor id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// if this instructor already exists
		if (c.instructors_map.find(str) != c.instructors_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - an instructor with this id was already defined"; return false; }
		
		// create the instructor
		instructor_t &inst = c.instructors.emplace_back();
		inst.id = str;
		c.instructors_map.emplace(std::move(str), &inst);

		if (!line_term() && f) { err = std::string(path) + ':' + std::to_string(ln) + " - unexpected tokens encountered after instructor id"; return false; }

		return true;
	}
	bool _parse_course()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected course id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// insert a new course into the set (make sure it doesn't already exist)
		if (c.courses_map.find(str) != c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to respecify existing course: " + str; return false; }
		course_t &course = c.courses.emplace_back();
		course.id = str;
		c.courses_map.emplace(std::move(str), &course);

		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected course capacity"; return false; }
		if (!(f >> course.capacity)) { err = std::string(path) + ':' + std::to_string(ln) + " - failed to parse course capacity"; return false; }

		// get the instructor for this course
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected instructor id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above
		auto inst = c.instructors_map.find(str);
		if (inst == c.instructors_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - reference to undefined instructor: " + str; return false; }
		course.instructor = inst->second; // fill in instructor field

		if (!line_term() && f) // if we have extra stuff it's notes for the course
		{
			std::getline(f, course.notes);
			assert(f); // guaranteed to succeed from above
			++line_term.line_number; // bump this up to account for getline() consuming the newline char
			if (auto p = course.notes.find('#'); p != std::string::npos) course.notes.erase(p); // remove comment from end of string (if present)
			course.notes.erase(course.notes.find_last_not_of(" \t\n\r\v\f") + 1); // trim white space from end of string
		}

		return true;
	}
	bool _parse_requires()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected course id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// find the course being constrained
		auto it = c.courses_map.find(str);
		if (it == c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str; return false; }

		// gather all the tokens
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected one or more room attributes"; return false; }
		do
		{
			f >> str;
			assert(f); // guaranteed to succeed from above
			it->second->required_attrs.insert(std::move(str)); // add requirement (duplicates are no-op)
		} while (!line_term() && f);

		return true;
	}
	bool _parse_parallel_orthogonal()
	{
		// we handle parallel and orthogonal with same code - use this to separate the semantics
		const auto bucket = str == "parallel" ? &course_t::parallel_courses : &course_t::orthogonal_courses;

		// parse the course list
		course_set.clear();
		while (!line_term() && f)
		{
			f >> str;
			assert(f); // guaranteed to succeed from above
			auto it = c.courses_map.find(str);
			if (it == c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str; return false; }
			course_set.insert(it->second);
		}
		if (course_set.size() < 2) { err = std::string(path) + ':' + std::to_string(ln) + " - expected two or more course ids"; return false; }

		// apply constraints
		for (auto *dest : course_set)
		{
			auto &b = dest->*bucket; // get the destination set
			for (auto *constraint : course_set) if (dest != constraint)
			{
				b.insert(constraint); // add all other constraints to dest
			}
		}

		return true;
	}
	bool _parse_follows()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected course id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// find the course being constrained (first)
		auto first = c.courses_map.find(str);
		if (first == c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str; return false; }

		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected a second course id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// find the course being constrained (second)
		auto second = c.courses_map.find(str);
		if (second == c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str; return false; }

		if (!line_term() && f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected only two arguments"; return false; }

		// apply the constraint
		first->second->follow_courses.insert(second->second);

		return true;
	}
	bool _parse_schedule()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected course id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// find the course being constrained
		auto it = c.courses_map.find(str);
		if (it == c.courses_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to constrain undefined course: " + str; return false; }

		// get the schedule id
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected a schedule id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above
		auto sched_ref = c.schedules_map.find(str);
		if (sched_ref == c.schedules_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - reference to undefined schedule: " + str; return false; }

		// if there's already a schedule constraint for this course and they differ, it's a problem
		if (schedule_t *old = it->second->required_schedule; old && old != sched_ref->second) { err = std::string(path) + ':' + std::to_string(ln) + " - attempt to override pre-existing schedule constraint: " + old->id + " -> " + str; return false; }

		// apply the constraint
		it->second->required_schedule = sched_ref->second;

		if (!line_term() && f) { err = std::string(path) + ':' + std::to_string(ln) + " - unexpected tokens encountered after schedule id"; return false; }

		return true;
	}
	bool _parse_unavailable()
	{
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected instructor id"; return false; }
		f >> str;
		assert(f); // guaranteed to succeed from above

		// get the instructor being modified - if not defined error
		auto inst = c.instructors_map.find(str);
		if (inst == c.instructors_map.end()) { err = std::string(path) + ':' + std::to_string(ln) + " - undefined reference to instructor: " + str; return false; }

		// parse the time
		timespan span;
		if (line_term() || !f) { err = std::string(path) + ':' + std::to_string(ln) + " - expected timespan after instructor id"; return false; }
		if (!parse_time_range(f, span)) { err = std::string(path) + ':' + std::to_string(ln) + " - failed to parse timespan"; return false; }
		if (!line_term() && f) { err = std::string(path) + ':' + std::to_string(ln) + " - unexpected tokens encountered after timespan"; return false; }

		// add it to unavailability
		inst->second->unavailability.add(span);

		return true;
	}
};
std::unordered_map<std::string, bool(constraint_parser_pack::*)()> parse_handlers
{
	{ "room", &constraint_parser_pack::_parse_room },
{ "timeslot", &constraint_parser_pack::_parse_timeslot },
{ "instructor", &constraint_parser_pack::_parse_instructor },
{ "course", &constraint_parser_pack::_parse_course },
{ "requires", &constraint_parser_pack::_parse_requires },
{ "parallel", &constraint_parser_pack::_parse_parallel_orthogonal },
{ "orthogonal", &constraint_parser_pack::_parse_parallel_orthogonal },
{ "follows", &constraint_parser_pack::_parse_follows },
{ "schedule", &constraint_parser_pack::_parse_schedule },
{ "unavailable", &constraint_parser_pack::_parse_unavailable },
};
std::variant<constraint_set, std::string> read_constraints(std::istream &f, const char *path)
{
	constraint_parser_pack p{ f, path };

	// parse all the course data
	while (true)
	{
		if (p.line_term()) continue; // if this is line term, line was empty (skip)
		if (!f) { assert(f.eof()); break; } // this happens at eof
		p.ln = p.line_term.line_number;

		f >> p.str;
		p.str = to_lower(std::move(p.str)); // convert instruction to lowercase to be case invariant
		assert(f); // guaranteed to succeed from above

		// get the parse handler - if it doesn't exist unknown command - otherwise execute
		auto it = parse_handlers.find(p.str);
		if (it == parse_handlers.end()) return std::string(path) + ':' + std::to_string(p.ln) + " - unrecognized command: " + p.str;
		if (!(p.*it->second)()) return std::move(p.err);
	}

	// ensure usage constraints for rooms
	if (p.c.rooms.empty()) return std::string(path) + " - no rooms were defined";

	// ensure usage constraints for schedules
	if (p.c.schedules.empty()) return std::string(path) + " - no schedules were defined";
	for (const auto &sched : p.c.schedules)
	{
		assert(!sched.slots.empty()); // this should be impossible

		// make sure no timeslots overlap
		for (std::size_t i = 1; i < sched.slots.size(); ++i)
		{
			if (sched.slots[i].time.start < sched.slots[i - 1].time.stop)
				return std::string(path) + " - schedule " + sched.id + " has overlapping time slots: " + tostr(sched.slots[i - 1].time) + " and " + tostr(sched.slots[i].time);
		}
	}

	return std::move(p.c);
}

struct satisfy_info
{
public:
	constraint_set &constraints;

	std::list<std::list<std::set<const course_t*>>>                              follow_chains;            // list of all follow chains in the course topology
	std::unordered_map<const std::list<std::set<const course_t*>>*, schedule_t*> follow_chain_to_schedule; // maps each follow chain to its required schedule or nullptr if no specific schedule is required

	// maps schedules to a map from follow chain to valid start index - used to severely prune state space for some types of constraints
	std::unordered_map<const schedule_t*, std::unordered_map<const std::list<std::set<const course_t*>>*, std::vector<std::size_t>>> schedule_to_follow_chain_to_start_index;

	std::unordered_map<const course_t*, std::vector<std::size_t>> course_to_satisfiable_rooms; // maps each course to the set of all its satisfiable rooms based on room attrs

	std::unordered_map<const std::set<const course_t*>*, std::set<const instructor_t*>> pset_to_instructors; // maps any pset to its set of instructors

	// constructs a new satisfiability info object
	// if this throws it means that the system is impossibly-constrained (as opposed to just being unsatisfiable)
	satisfy_info(constraint_set &_c) : constraints(_c)
	{
		std::list<std::set<const course_t*>>                                                   psets;                  // list of all parallel course sets
		std::unordered_map<const course_t*, std::list<std::set<const course_t*>>::iterator>    course_to_pset;         // maps each course to its pset
		std::unordered_map<const std::set<const course_t*>*, const std::set<const course_t*>*> follow_psets;           // maps a pset to its follow pset (if any)
		std::vector<std::list<std::set<const course_t*>>::iterator>                            fset;                   // temporary for building follow sets
		bool                                                                                   fmerge;                 // flag for looping construction logic
		std::unordered_map<const std::set<const course_t*>*, time_union>                       pset_to_unavailability; // maps each pset to its total unavailability union

		// construct the trivial form of the parallel topology structure - just a bunch of bookmarked singletons
		for (const auto &entry : constraints.courses)
		{
			psets.emplace_back().insert(&entry);
			course_to_pset.emplace(&entry, std::prev(psets.end()));
		}

		// now apply all parallel constraints
		for (auto pos = psets.begin(); pos != psets.end(); )
		{
			bool merged = false;

			// for each course in the current parallel set
			for (const auto *p : *pos)
			{
				// for each parallel constraint
				for (const auto *pp : p->parallel_courses)
				{
					// get its parallel set - if it's us skip it
					auto it = course_to_pset.at(pp);
					if (it == pos) continue;

					// merge it into our parallel set and account for all mapping updates
					for (const auto *ppp : *it) course_to_pset.at(ppp) = pos;
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
					for (const auto *p : *other)
					{
						for (const auto *pp : p->follow_courses)
						{
							// get the follow pset
							auto val = course_to_pset.at(pp);
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
					for (const auto *g : *fset[i]) course_to_pset.at(g) = dest;
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
				for (const auto *i : pset)
				{
					for (const auto *j : i->follow_courses)
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
					for (const auto *g : *fset[i]) course_to_pset.at(g) = dest;
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
			schedule_t     *required_schedule = nullptr; // keep track of the required schedule for this follow chain
			const course_t *source = nullptr;            // also keep track of the source for this constraint for error messages
			
			for (const auto &pset : chain)
			{
				for (const auto *p : pset)
				{
					// get this course - if this course doesn't require a specific schedule skip it
					if (!p->required_schedule) continue;

					// if there's a required schedule for this chain we need to match it
					if (required_schedule)
					{
						if (required_schedule != p->required_schedule) throw std::logic_error(std::string("conflicting schedule specifiers in same follow chain: ") + source->id + " -> " + required_schedule->id + " vs " + p->id + " -> " + p->required_schedule->id);
					}
					// otherwise just mark this as the new required schedule for the chain
					else
					{
						required_schedule = p->required_schedule;
						source = p;
					}
				}
			}

			follow_chain_to_schedule.emplace(&chain, required_schedule); // add this chain to the map
		}

		// now we generate the total unavailability unions for each pset, as well as the instructor sets
		for (const auto &chain : follow_chains)
		{
			for (const auto &pset : chain)
			{
				time_union                    unavail;
				std::set<const instructor_t*> instructors;

				for (const auto *i : pset)
				{
					unavail.add(i->instructor->unavailability);
					instructors.insert(i->instructor);
				}

				pset_to_unavailability.emplace(&pset, std::move(unavail));
				pset_to_instructors.emplace(&pset, std::move(instructors));
			}
		}

		// examine each schedule
		for (const auto &sched : constraints.schedules)
		{
			std::unordered_map<const std::list<std::set<const course_t*>>*, std::vector<std::size_t>> follow_chain_to_start_index;

			// and each follow chain
			for (const auto &chain : follow_chains)
			{
				std::vector<std::size_t> good_starts;

				// and each starting position
				for (std::size_t start = 0; start < sched.slots.size() && start + chain.size() <= sched.slots.size(); ++start)
				{
					bool good = true;

					// and each pset in the chain
					auto pset = chain.begin();
					for (std::size_t i = 0; i < chain.size(); ++i, ++pset)
					{
						// if this violates an availability, no good
						if (pset_to_unavailability.at(&*pset).intersects(sched.slots[start + i].time))
						{
							good = false;
							break;
						}
					}

					// if this starting position was good, insert it as a valid starting 
					if (good) good_starts.push_back(start);
				}

				// if there were no good starts this is impossible
				if (good_starts.empty()) throw std::logic_error("no valid starting positions for one or more follow chains");
				// otherwise insert into the follow chain map
				follow_chain_to_start_index.emplace(&chain, std::move(good_starts));
			}

			// insert into the schedule map
			schedule_to_follow_chain_to_start_index.emplace(&sched, std::move(follow_chain_to_start_index));
		}

		// build map from each course to its set of satisfiable rooms
		for (const auto &course : constraints.courses)
		{
			std::vector<std::size_t> valid_rooms;

			for (std::size_t i = 0; i < constraints.rooms.size(); ++i)
			{
				if (constraints.rooms[i].capacity < course.capacity) continue;

				const auto &attrs = constraints.rooms[i].attr;

				bool valid = true;
				for (const auto &attr : course.required_attrs)
				{
					if (attrs.find(attr) == attrs.end())
					{
						valid = false;
						break;
					}
				}
				if (valid) valid_rooms.push_back(i);
			}

			if (valid_rooms.empty()) throw std::logic_error(std::string("course ") + course.id + " has no valid rooms due to required attr constraints");

			course_to_satisfiable_rooms.emplace(&course, std::move(valid_rooms));
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
				ostr << "    ";
				for (const auto *ss : s) ostr << ' ' << std::setw(16) << ss->id;
				ostr << '\n';
			}
		}
		ostr << '\n';
	}

	// begins by erasing the current schedules, then recursively applies constraints in order to generate a valid schedule.
	// returns true on success (and leaves the result in schedule).
	bool satisfy();

private:
	bool satisfy_fchain(std::list<std::list<std::set<const course_t*>>>::const_iterator chain);
	bool satisfy_pset_at(std::list<std::list<std::set<const course_t*>>>::const_iterator chain, schedule_t &sched, std::size_t slot_i, std::list<std::set<const course_t*>>::const_iterator pset);
	bool satisfy_pset_at_interior(std::list<std::list<std::set<const course_t*>>>::const_iterator chain, schedule_t &sched, std::size_t slot_i, std::list<std::set<const course_t*>>::const_iterator pset, std::set<const course_t*>::const_iterator ppos);
};

bool satisfy_info::satisfy_pset_at_interior(std::list<std::list<std::set<const course_t*>>>::const_iterator chain, schedule_t &sched, std::size_t slot_i, std::list<std::set<const course_t*>>::const_iterator pset, std::set<const course_t*>::const_iterator ppos)
{
	// if we're at the end of the current pset, we're done with this pset
	if (ppos == pset->end())
	{
		// recurse to the next item in the follow chain
		return satisfy_pset_at(chain, sched, slot_i + 1, std::next(pset));
	}

	const auto   &course = **ppos;
	timeslot_t   &slot = sched.slots[slot_i];

	// attempt to put it into each available room
	assert(slot.assignments.size() == constraints.rooms.size());
	for (std::size_t k : course_to_satisfiable_rooms.at(&course))
	{
		// if this room is already taken, it's not viable
		if (slot.assignments[k]) continue;

		// if we violate an orthogonality constraint, it's not viable
		if ([&] {
			const auto &ortho = course.orthogonal_courses;
				for (const auto *other_assignment : slot.assignments)
				{
					if (ortho.find(const_cast<course_t*>(other_assignment)) != ortho.end()) return true;
				}
			return false;
		}()) continue;

		// ----------------------------------------------------------------

		// assign current course to this schedule, timeslot, and room
		slot.assignments[k] = *ppos;

		// perform the recursive step - if we succeed, propagate up
		if (satisfy_pset_at_interior(chain, sched, slot_i, pset, std::next(ppos))) return true;

		// otherwise undo the change and continue searching
		slot.assignments[k] = nullptr;
	}

	return false;
}
bool satisfy_info::satisfy_pset_at(std::list<std::list<std::set<const course_t*>>>::const_iterator chain, schedule_t &sched, std::size_t slot_i, std::list<std::set<const course_t*>>::const_iterator pset)
{
	// if we're at the end of the follow chain, we're done with this chain - recurse to the next
	if (pset == chain->end()) return satisfy_fchain(std::next(chain));

	assert(slot_i < sched.slots.size()); // sanity check
	timeslot_t &slot = sched.slots[slot_i];
	const auto &instructors = pset_to_instructors.at(&*pset);

	// if this timeslot is already scheduled for any instructor in this pset it's no good
	for (const auto *ins : slot.assignments)
	{
		if (ins && instructors.find(ins->instructor) != instructors.end()) return false;
	}

	// otherwise recurse into the pset
	return satisfy_pset_at_interior(chain, sched, slot_i, pset, pset->begin());
}
bool satisfy_info::satisfy_fchain(std::list<std::list<std::set<const course_t*>>>::const_iterator chain)
{
	// if we're at the end of the follow chains we're done and have scheduled everything successfully (yay)
	if (chain == follow_chains.end()) return true;

	// create a function to attempt to place the current course into the given schedule
	auto scheduler = [&] (schedule_t &sched) {
		// try to put it into each valid starting timeslot in the given schedule
		for (std::size_t slot_i : schedule_to_follow_chain_to_start_index.at(&sched).at(&*chain))
		{
			// if there's not enough room to fit the entire chain we can prune this entire execution branch
			//if (slot_i + chain->size() > sched.slots.size()) break;
			assert(slot_i + chain->size() <= sched.slots.size());

			// if we can satisfy this pset with this slot
			if (satisfy_pset_at(chain, sched, slot_i, chain->begin())) return true;
		}
		return false;
	};

	// get the required schedule for this chain - if there is one, use that
	if (schedule_t *req = follow_chain_to_schedule.at(&*chain); req)
	{
		return scheduler(*req);
	}
	// otherwise attempt to put current follow chain into each schedule in order
	else
	{
		for (auto &sched : constraints.schedules)
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
	for (auto &sched : constraints.schedules)
	{
		for (auto &slot : sched.slots)
		{
			slot.assignments.clear();
			slot.assignments.resize(constraints.rooms.size());
		}
	}

	// recurse to satisfy all course assignments
	return satisfy_fchain(follow_chains.begin());
}

void print_schedule_latex(std::ostream &ostr, const constraint_set &c, const schedule_t &sched)
{
	ostr << "\\begin{table}[ht!]\n\\centering\n\\begin{tabularx}{\\textwidth}{|X";
	for (std::size_t i = 0; i < c.rooms.size(); ++i) ostr << "|X";
	ostr << "|}\n\\hline ";

	ostr << fix_id(sched.id);
	for (const auto &room : c.rooms) ostr << " & " << fix_id(room.id);
	ostr << " \\\\\n\\hline ";

	for (const auto &timeslot : sched.slots)
	{
		ostr << timeslot.time;
		for (const auto *asgn : timeslot.assignments)
		{
			ostr << " & ";
			if (asgn)
			{
				ostr << fix_id(asgn->id) << ' ' << fix_id(asgn->instructor->id);
				if (!asgn->notes.empty()) ostr << ' ' << asgn->notes;
			}
		}
		ostr << " \\\\\n\\hline ";
	}

	ostr << "\n\\end{tabularx}\n\\end{table}\n";
}
void print_schedules_latex(std::ostream &ostr, const constraint_set &c)
{
	ostr << "\\documentclass[8pt]{article}\n\\usepackage[margin=0.5in]{geometry}\n\\usepackage{tabularx}\n\\begin{document}\n\n";

	for (const auto &sched : c.schedules)
	{
		print_schedule_latex(ostr, c, sched);
		ostr << '\n';
	}

	ostr << "\\end{document}\n";
}

void print_schedule_text(std::ostream &ostr, const constraint_set &c, const schedule_t &sched)
{
	ostr << fix_id(sched.id);
	for (const auto &room : c.rooms) ostr << '\t' << fix_id(room.id);
	ostr << '\n';

	for (const auto &timeslot : sched.slots)
	{
		ostr << timeslot.time;
		for (const auto *asgn : timeslot.assignments)
		{
			ostr << '\t';
			if (asgn)
			{
				ostr << fix_id(asgn->id) << ' ' << fix_id(asgn->instructor->id);
				if (!asgn->notes.empty()) ostr << ' ' << asgn->notes;
			}
		}
		ostr << '\n';
	}
}
void print_schedules_text(std::ostream &ostr, const constraint_set &c)
{
	for (const auto &sched : c.schedules)
	{
		print_schedule_text(ostr, c, sched);
		ostr << '\n';
	}
}

const char *const help_msg = R"(Usage: schedule [OPTION]... [constraints file]
Generate class schedule given rooms, timeslots, and course info.
If constraints file is not specified reads from stdin.

  --help               print this help page and exit
  --latex              generate latex source instead of tab-separated text
  --topology           print the course topology to stderr in addition to all other work
)";

int main(int argc, const char *const argv[]) try
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

	if (pathspec.size() > 1)
	{
		std::cerr << "incorrect usage: see --help for info\n";
		return 1;
	}

	// ---------------------------------------------

	// parse the input file (or stdin if none was specified)
	std::variant<constraint_set, std::string> _constraints;
	if (!pathspec.empty())
	{
		std::ifstream f{ pathspec[0] };
		if (!f)
		{
			std::cerr << "failed to open " << pathspec[0] << " for reading\n";
			return 60;
		}
		_constraints = read_constraints(f, pathspec[0]);
	}
	else _constraints = read_constraints(std::cin, "<stdin>");

	// ---------------------------------------------

	// check the error variant for parse result
	if (_constraints.index() != 0)
	{
		std::cerr << std::get<1>(_constraints) << '\n';
		return 100;
	}
	auto &constraints = std::get<0>(_constraints);

	// ---------------------------------------------

	// attempt to create the solver object - if this fails it means that there were impossible constraints
	std::unique_ptr<satisfy_info> solver;
	try { solver = std::make_unique<satisfy_info>(constraints); }
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
	printer(std::cout, constraints);

	return 0;
}
catch (const std::exception & e)
{
	std::cerr << "ERROR: unhandled exception\n" << e.what() << '\n';
	return 600;
}
catch (...)
{
	std::cerr << "ERROR: unhandled exception of unknown type\n";
	return 601;
}
