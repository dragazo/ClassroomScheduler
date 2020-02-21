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
		else return std::string(path) + ':' + std::to_string(ln) + " - unrecognized command: " + str;
	}

	return std::move(courses);
}

struct satisfy_info
{
	const std::vector<room> &rooms;
	std::vector<schedule> &schedules;
	const std::unordered_map<std::string, course_info> &courses;

	std::list<std::set<std::string>> parallel_sets; // list of parallel course sets

	satisfy_info(const std::vector<room> &r, std::vector<schedule> &s, const std::unordered_map<std::string, course_info> &c) : rooms(r), schedules(s), courses(c)
	{
		// construct the parallel topology object for generative iteration
		for (const auto &entry : c)
		{
			auto p_pos = entry.second.parallel_courses.begin();
			const auto p_end = entry.second.parallel_courses.end();
			bool inserted = false;

			// for the current course, iterate through all its parallel courses (p)
			for (; p_pos != p_end; ++p_pos)
			{
				// look for a parallel set that contains p
				auto it = std::find_if(parallel_sets.begin(), parallel_sets.end(), [&](const auto &x) { return x.find(*p_pos) != x.end(); });
				if (it != parallel_sets.end())
				{
					// insert ourself into the matching parallel set
					it->insert(entry.first);
					inserted = true;

					// continue iterating through all the rest of the parallel courses (p) - we know previous ones failed to match any parallel sets
					for (++p_pos; p_pos != p_end; ++p_pos)
					{
						// look for anoter parallel set that contains p
						auto other = std::find_if(parallel_sets.begin(), parallel_sets.end(), [&](const auto &x) { return x.find(*p_pos) != x.end(); });
						if (other != parallel_sets.end() && other != it)
						{
							// if we found another match elsewhere, merge into this one and erase it - this gives us transitive parallel constraints
							it->merge(*other);
							parallel_sets.erase(other);
							std::cout << "performed a merge!\n";
						}
					}

					break;
				}
			}
			// if we reached the end of the container then there was no matching parallel set - create one and put ourself in it
			if (!inserted) parallel_sets.emplace_back().insert(entry.first);
		}

		//parallel_sets.sort([](const auto &a, const auto &b) { return a.size() > b.size(); });
	}
};
bool _satisfy_recursive(satisfy_info &p, std::list<std::set<std::string>>::const_iterator pset, timeslot &slot, std::set<std::string>::const_iterator pset_item)
{
	// if we're at the end of the current pset, we're done with this pset
	assert(pset != p.parallel_sets.end());
	if (pset_item == pset->end()) return true;

	const auto &course = p.courses.at(*pset_item);

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
		slot.assignments[k] = *pset_item;

		// perform the recursive step - if we succeed, propagate up
		if (_satisfy_recursive(p, pset, slot, std::next(pset_item))) return true;

		// otherwise undo the change and continue searching
		slot.assignments[k].clear();
	}

	return false;
}
bool _satisfy_recursive(satisfy_info &p, std::list<std::set<std::string>>::const_iterator pset)
{
	// if we're at the end of the parallel sets, we're done
	if (pset == p.parallel_sets.end()) return true;

	// attempt to put current course set into each schedule
	for (auto &sched : p.schedules)
	{
		// and into each timeslot in said schedule
		for (auto &slot : sched.slots)
		{
			// if we can setisfy this pset with this slot
			if (_satisfy_recursive(p, pset, slot, pset->begin()))
			{
				// recurse to handle the next pset - if it succeeds, propagate back up
				if (_satisfy_recursive(p, std::next(pset))) return true;
			}
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

	// create info object and begin recursive solve strategy
	satisfy_info p{ rooms, schedules, courses };
	return _satisfy_recursive(p, p.parallel_sets.begin());
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